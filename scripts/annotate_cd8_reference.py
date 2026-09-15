"""Sample-independent external-reference CD8 annotation; preserves v1 analyses."""
import argparse
import gc
import gzip
import hashlib
import json
import os
import shutil
import time
from pathlib import Path
os.environ.setdefault('OMP_NUM_THREADS','4')
os.environ.setdefault('NUMBA_NUM_THREADS','4')
_ROOT=Path(__file__).resolve().parents[1]
_CACHE=_ROOT/'models/celltypist_cache'
(_CACHE/'data/models').mkdir(parents=True,exist_ok=True)
for _name in ['Immune_All_Low.pkl','Immune_All_High.pkl']:
    if not (_CACHE/'data/models'/_name).exists():
        shutil.copy2(_ROOT/'models/celltypist'/_name,_CACHE/'data/models'/_name)
os.environ['CELLTYPIST_FOLDER']=str(_CACHE)
import anndata as ad
import celltypist
from celltypist.models import Model
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmread
from threadpoolctl import threadpool_limits
from run_decoded_state_baselines import PANELS, RAW, ROOT

OUT=ROOT/'results/cd8_state_v2'
OBJECTS=ROOT/'data/processed/cd8_state_v2'
CD8_LABELS={'Tcm/Naive cytotoxic T cells','Tem/Temra cytotoxic T cells',
            'Tem/Trm cytotoxic T cells','Trm cytotoxic T cells','CD8a/b(entry)'}
MARKERS=['CD3D','CD3E','TRAC','CD8A','CD8B','CD4','FOXP3','NKG7','GNLY',
         'TRDC','TRGC1','CD79A','MS4A1','MZB1','LYZ','LST1','FCN1','S100A8',
         'EPCAM','KRT18','KRT19']

def annotate_sample(patient,point,low,high):
    sid=f'{patient}_{point}'
    sample_dir=OUT/sid
    sample_dir.mkdir(parents=True,exist_ok=True)
    done=sample_dir/'summary.json'
    if done.exists():
        return json.loads(done.read_text(encoding='utf-8'))
    start=time.perf_counter()
    stem=RAW/f'GSE288199_{sid}_T'
    with gzip.open(str(stem)+'_matrix.mtx.gz','rb') as f:
        counts=mmread(f,spmatrix=True).T.tocsr().astype(np.float32)
    counts.sum_duplicates(); counts.eliminate_zeros()
    with gzip.open(str(stem)+'_features.tsv.gz','rt') as f:
        features=[l.rstrip().split('\t') for l in f]
    with gzip.open(str(stem)+'_barcodes.tsv.gz','rt') as f:
        barcodes=[l.strip() for l in f]
    assert counts.shape==(len(barcodes),len(features))
    assert len(barcodes)==len(set(barcodes))
    genes=np.array([r[1] for r in features])
    total=np.asarray(counts.sum(axis=1)).ravel()
    ngenes=np.diff(counts.indptr)
    mt=np.asarray(counts[:,np.char.startswith(genes,'MT-')].sum(axis=1)).ravel()
    mtpct=np.divide(mt*100,total,out=np.zeros_like(total),where=total>0)
    keep=(total>=500)&(ngenes>=200)&(mtpct<25)
    unique,inverse=np.unique(genes,return_inverse=True)
    if len(unique)!=len(genes):
        collapse=sparse.csr_matrix((np.ones(len(genes)),(np.arange(len(genes)),inverse)),
                                  shape=(len(genes),len(unique)),dtype=np.float32)
        counts=(counts@collapse).tocsr()
        genes=unique
    counts=counts[keep].copy()
    obs=pd.DataFrame({'barcode':np.array(barcodes)[keep], 'patient_id':patient,'timepoint':point,
                      'total_counts':total[keep],'n_genes':ngenes[keep],'pct_mito':mtpct[keep]},
                     index=[f'{sid}_{b}' for b in np.array(barcodes)[keep]])
    raw=ad.AnnData(counts,obs=obs,var=pd.DataFrame(index=genes))
    print(f'{sid}: QC {raw.n_obs}, doublet screen',flush=True)
    sc.pp.scrublet(raw,expected_doublet_rate=.05,random_state=20260909,n_prin_comps=30,
                   use_approx_neighbors=False,verbose=False)
    assert 'predicted_doublet' in raw.obs and raw.obs.predicted_doublet.notna().all()
    assert np.isfinite(raw.obs.doublet_score).all()
    norm=raw.copy()
    sc.pp.normalize_total(norm,target_sum=10000)
    sc.pp.log1p(norm)
    labels=[]
    for first in range(0,norm.n_obs,2048):
        block=norm[first:first+2048].copy()
        lo=celltypist.annotate(block,model=low,majority_voting=False)
        hi=celltypist.annotate(block,model=high,majority_voting=False)
        labels.append(pd.DataFrame({'reference_low':lo.predicted_labels['predicted_labels'].astype(str),
                                    'reference_confidence':lo.probability_matrix.max(axis=1),
                                    'reference_high':hi.predicted_labels['predicted_labels'].astype(str)}))
    labels=pd.concat(labels).loc[raw.obs_names]
    for col in labels:
        raw.obs[col]=labels[col]
    present=lambda panel: np.asarray((counts[:,np.isin(genes,panel)]>0).sum(axis=1)).ravel()
    marker_t=present(['CD3D','CD3E','TRAC'])>0
    marker_cd8=present(['CD8A','CD8B'])>0
    raw.obs['old_cd8_like']=marker_t&(present(['CD8B'])>0)
    raw.obs['mixed_lineage_flag']=(present(['CD79A','MS4A1','MZB1'])>=2)|(
        present(['LYZ','LST1','FCN1','S100A8'])>=2)|(present(['EPCAM','KRT18','KRT19'])>=2)
    selected=(raw.obs.reference_low.isin(CD8_LABELS).to_numpy() &
              (raw.obs.reference_confidence.to_numpy()>=.5)&
              (raw.obs.reference_high.to_numpy()=='T cells')&marker_t&marker_cd8&
              ~raw.obs.predicted_doublet.to_numpy())
    raw.obs['reference_cd8_selected']=selected
    for g in MARKERS:
        raw.obs['expr_'+g]=np.asarray(norm.X[:,genes==g].sum(axis=1)).ravel()
    raw.obs.to_csv(sample_dir/'cell_annotations.tsv.gz',sep='\t')
    raw.obs.groupby('reference_low',observed=True).agg(n=('barcode','size'),
        selected=('reference_cd8_selected','sum'),confidence=('reference_confidence','mean')).to_csv(
        sample_dir/'reference_counts.tsv',sep='\t')
    marker_table=raw.obs.groupby('reference_low',observed=True)[['expr_'+g for g in MARKERS]].mean()
    marker_table.to_csv(sample_dir/'marker_means.tsv',sep='\t')
    score_rows=[]
    for compartment,mask in [('all_qc',np.ones(raw.n_obs,dtype=bool)),('cd8_like',selected)]:
        row=dict(patient_id=patient,timepoint=point,compartment=compartment,
                 raw_cells=len(barcodes),qualifying_cells=int(mask.sum()),eligible=int(mask.sum())>=50)
        for panel,markers in PANELS.items():
            assert set(markers).issubset(set(genes))
            row[panel]=float(norm.X[mask][:,np.isin(genes,markers)].mean()) if mask.any() else np.nan
        score_rows.append(row)
    pd.DataFrame(score_rows).to_csv(sample_dir/'program_scores.tsv',sep='\t',index=False)
    # Persist the complete selected pretreatment population for reproducible sampling/embedding.
    if point=='pre' and selected.any():
        chosen=raw[selected].copy()
        chosen.var['gene_symbols']=chosen.var_names.astype(str)
        chosen.write_h5ad(OBJECTS/f'{sid}_reference_cd8_counts.h5ad',compression='gzip')
    summary=dict(patient_id=patient,timepoint=point,raw_cells=len(barcodes),qc_cells=raw.n_obs,
                 doublets=int(raw.obs.predicted_doublet.sum()),
                 scrublet_threshold=float(raw.uns['scrublet']['threshold']),
                 reference_cd8=int(selected.sum()),old_cd8_like=int(raw.obs.old_cd8_like.sum()),
                 selected_mixed_lineage=int(raw.obs.loc[selected,'mixed_lineage_flag'].sum()),
                 old_selected_overlap=int((selected&raw.obs.old_cd8_like.to_numpy()).sum()),
                 seconds=round(time.perf_counter()-start,2))
    done.write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary),flush=True)
    del raw,norm,counts
    gc.collect()
    return summary

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--patients',nargs='*')
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True); OBJECTS.mkdir(parents=True,exist_ok=True)
    low=Model.load(str(ROOT/'models/celltypist/Immune_All_Low.pkl'))
    high=Model.load(str(ROOT/'models/celltypist/Immune_All_High.pkl'))
    assert CD8_LABELS.issubset(set(low.cell_types))
    manifest=pd.read_csv(ROOT/'results/decoded_state_pilot/patient_manifest.tsv',sep='\t')
    patients=manifest.loc[manifest.candidate_pair,'patient_id'].tolist()
    if args.patients:
        assert set(args.patients).issubset(set(patients))
        patients=args.patients
    summaries=[]
    with threadpool_limits(limits=4):
        for patient in patients:
            for point in ['pre','post']:
                summaries.append(annotate_sample(patient,point,low,high))
    pd.DataFrame(summaries).to_csv(OUT/'annotation_summary.tsv',sep='\t',index=False)
    print('Annotation finished',flush=True)

if __name__=='__main__':
    main()

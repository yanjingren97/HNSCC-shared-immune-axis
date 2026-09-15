"""Patient-level composition audit and summed raw counts; no outcomes loaded."""
from pathlib import Path
import gzip, json, hashlib, time
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/signature_revision_20260915/composition'
CFG_PATH=ROOT/'config/signature_composition_revision_20260915.json'
CFG=json.loads(CFG_PATH.read_text(encoding='utf-8'))
GROUPS=CFG['primary_groups']
CD8=set(CFG['cd8_labels'])
T_LABELS=set().union(*map(set,GROUPS.values()),CD8,set(CFG['additional_T_labels']))

def annotation(patient):
    path=ROOT/f'results/cd8_state_v2/{patient}_pre/cell_annotations.tsv.gz'
    if not path.exists():
        path=ROOT/f'results/signature_revision_20260915/additional_annotations/{patient}_pre/cell_annotations.tsv.gz'
    return path

def masks(obs):
    def flag(c): return obs[c].astype(str).str.lower().eq('true').to_numpy()
    marker_t=obs[['expr_CD3D','expr_CD3E','expr_TRAC']].max(axis=1).gt(0).to_numpy()
    marker_cd8=obs[['expr_CD8A','expr_CD8B']].max(axis=1).gt(0).to_numpy()
    base=obs.reference_confidence.ge(CFG['confidence_min']).to_numpy() & obs.reference_high.eq('T cells').to_numpy() & marker_t & ~flag('predicted_doublet')
    concordant=base & obs.reference_low.isin(T_LABELS).to_numpy()
    clean=concordant & ~flag('mixed_lineage_flag')
    out={'all_QC':np.ones(len(obs),bool),'T_including_mixed':concordant,'T_clean':clean,
         'CD8_clean':clean & obs.reference_low.isin(CD8).to_numpy() & marker_cd8}
    for name,labels in GROUPS.items():
        out[name]=clean & obs.reference_low.isin(labels).to_numpy()
        if name.startswith('CD8'): out[name]&=marker_cd8
    return out,dict(label_conflicts=int((base&~obs.reference_low.isin(T_LABELS).to_numpy()).sum()),
        mixed_T=int((concordant&flag('mixed_lineage_flag')).sum()))

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'counts').mkdir(exist_ok=True)
    # Patient identity list only: expression matrix header contains no response values.
    patients=pd.read_csv(ROOT/'results/published_signature_external_audit/main_full_pretreatment_T_pseudobulk_logcpm.tsv.gz',sep='\t',nrows=0).columns[1:].tolist()
    qcrows=[]; countrows=[]; sources=[]
    for patient in sorted(patients):
        start=time.time(); dest=OUT/f'counts/{patient}.npz'
        ann=annotation(patient); obs=pd.read_csv(ann,sep='\t',index_col=0)
        assert obs.barcode.is_unique and obs.patient_id.eq(patient).all() and obs.timepoint.eq('pre').all()
        stem=ROOT/f'data/raw/GSE288199_decoded/GSE288199_{patient}_pre_T'
        barcodes=pd.read_csv(str(stem)+'_barcodes.tsv.gz',sep='\t',header=None)[0].astype(str)
        assert barcodes.is_unique
        idx=pd.Index(barcodes).get_indexer(obs.barcode)
        assert (idx>=0).all() and len(set(idx))==len(idx)
        select,extra=masks(obs)
        if not dest.exists():
            with gzip.open(str(stem)+'_matrix.mtx.gz','rb') as f: raw=mmread(f).tocsr()
            raw.sum_duplicates(); raw.eliminate_zeros()
            features=pd.read_csv(str(stem)+'_features.tsv.gz',sep='\t',header=None)
            genes=features.iloc[:,1].astype(str).str.upper().to_numpy()
            assert raw.shape==(len(genes),len(barcodes))
            assert np.allclose(np.asarray(raw[:,idx].sum(axis=0)).ravel(),obs.total_counts.to_numpy())
            names=list(select)
            values=np.column_stack([np.asarray(raw[:,idx[m]].sum(axis=1)).ravel() for m in select.values()])
            unique,inverse=np.unique(genes,return_inverse=True)
            collapse=sparse.csr_matrix((np.ones(len(genes)),(inverse,np.arange(len(genes)))),shape=(len(unique),len(genes)))
            values=np.asarray(collapse@values,dtype=np.float64)
            np.savez_compressed(dest,genes=unique.astype(str),groups=np.asarray(names),counts=values,
                n_cells=np.array([m.sum() for m in select.values()]),
                config_sha256=hashlib.sha256(CFG_PATH.read_bytes()).hexdigest(),
                annotation_sha256=hashlib.sha256(ann.read_bytes()).hexdigest())
            del raw,values
        cache=np.load(dest)
        assert str(cache['config_sha256'])==hashlib.sha256(CFG_PATH.read_bytes()).hexdigest()
        assert str(cache['annotation_sha256'])==hashlib.sha256(ann.read_bytes()).hexdigest()
        for j,name in enumerate(cache['groups']):
            assert int(cache['n_cells'][j])==int(select[name].sum())
            countrows.append(dict(patient_id=patient,compartment=name,n_cells=int(cache['n_cells'][j]),
                total_UMI=float(cache['counts'][:,j].sum()),fraction_of_clean_T=float(select[name].sum()/max(1,select['T_clean'].sum()))))
        qcrows.append(dict(patient_id=patient,raw_cells=len(barcodes),annotated_QC_cells=len(obs),
            clean_T=int(select['T_clean'].sum()),clean_CD8=int(select['CD8_clean'].sum()),**extra))
        sources.append(dict(patient_id=patient,annotation=str(ann.relative_to(ROOT)),annotation_sha256=str(cache['annotation_sha256']),barcode_alignment='unique exact match; count sums checked'))
        print(f'{patient}: T={select["T_clean"].sum()}, CD8={select["CD8_clean"].sum()}, {time.time()-start:.1f}s',flush=True)
    pd.DataFrame(qcrows).to_csv(OUT/'patient_QC.tsv',sep='\t',index=False)
    pd.DataFrame(countrows).to_csv(OUT/'compartment_counts.tsv',sep='\t',index=False)
    pd.DataFrame(sources).to_csv(OUT/'annotation_sources.tsv',sep='\t',index=False)
    counts=pd.DataFrame(countrows).pivot(index='patient_id',columns='compartment',values='n_cells')
    audit=dict(n_baseline_patients=len(patients),matched_four_groups={str(k):counts.index[(counts[list(GROUPS)]>=k).all(axis=1)].tolist() for k in CFG['minimum_cells']},
        config_sha256=hashlib.sha256(CFG_PATH.read_bytes()).hexdigest(),outcomes_used=False)
    (OUT/'build_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(json.dumps(audit,indent=2))

if __name__=='__main__':
    with threadpool_limits(limits=4): main()

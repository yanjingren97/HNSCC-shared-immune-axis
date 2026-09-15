"""Download verified deposited triplets and run the frozen exploratory baseline screen."""
import csv
import gzip
import hashlib
import json
import shutil
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/decoded_state_pilot'
RAW = ROOT/'data/raw/GSE288199_decoded'
PANELS = {
    'cytotoxic':['PRF1','GZMB','GNLY','NKG7','CTSW'],
    'interferon':['ISG15','IFIT1','IFIT2','IFIT3','MX1','OAS1'],
    'exhaustion_associated':['PDCD1','LAG3','HAVCR2','TIGIT','TOX','CTLA4'],
}

def fetch(row):
    name = row['url'].rsplit('/',1)[1]
    target = RAW/name
    existing = ROOT/'data/pilot/GSE288199_HN01'/name
    if not target.exists() and existing.exists():
        shutil.copy2(existing,target)
    if not target.exists():
        if shutil.disk_usage(RAW).free < 5*1024**3:
            raise RuntimeError('Less than 5 GiB free; download stopped')
        partial = target.with_suffix(target.suffix+'.part')
        with urllib.request.urlopen(row['url'],timeout=60) as response, partial.open('wb') as f:
            shutil.copyfileobj(response,f)
        # Full gzip CRC validation before promoting the file.
        with gzip.open(partial,'rb') as f:
            while f.read(1024*1024):
                pass
        partial.replace(target)
    return dict(**row,bytes=target.stat().st_size,sha256=hashlib.sha256(target.read_bytes()).hexdigest())

def score(patient,timepoint):
    stem=RAW/f'GSE288199_{patient}_{timepoint}_T'
    with gzip.open(str(stem)+'_features.tsv.gz','rt') as f:
        features=[line.rstrip().split('\t') for line in f]
    with gzip.open(str(stem)+'_barcodes.tsv.gz','rt') as f:
        barcodes=[line.strip() for line in f]
    with gzip.open(str(stem)+'_matrix.mtx.gz','rb') as f:
        matrix=mmread(f).T.tocsr().astype(np.float32)
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    assert matrix.shape==(len(barcodes),len(features))
    assert len(set(barcodes))==len(barcodes)
    assert np.isfinite(matrix.data).all() and (matrix.data>=0).all()
    assert np.equal(matrix.data,np.floor(matrix.data)).all()
    genes=np.array([r[1] for r in features])
    required=set(sum(PANELS.values(),[])+['CD8B','CD3D','CD3E','TRAC'])
    assert required.issubset(set(genes)), required-set(genes)
    assert all(sum(genes==g)==1 for g in required), 'Ambiguous marker gene symbols'
    total=np.asarray(matrix.sum(axis=1)).ravel()
    detected=np.diff(matrix.indptr)
    mito=np.asarray(matrix[:,np.char.startswith(genes,'MT-')].sum(axis=1)).ravel()
    pct=np.divide(100*mito,total,out=np.zeros_like(total),where=total>0)
    keep=(total>=500)&(detected>=200)&(pct<25)
    cd8=(np.asarray(matrix[:,genes=='CD8B'].sum(axis=1)).ravel()>0)
    tcell=(np.asarray(matrix[:,np.isin(genes,['CD3D','CD3E','TRAC'])].sum(axis=1)).ravel()>0)
    masks={'all_qc':keep,'cd8_like':keep&cd8&tcell}
    norm=matrix.multiply(np.divide(10000,total,out=np.zeros_like(total),where=total>0)[:,None]).tocsr()
    norm.data=np.log1p(norm.data)
    rows=[]
    for compartment,mask in masks.items():
        r=dict(patient_id=patient,timepoint=timepoint,compartment=compartment,
               raw_cells=len(barcodes),qualifying_cells=int(mask.sum()),eligible=int(mask.sum())>=50)
        for panel,markers in PANELS.items():
            r[panel]=float(norm[mask][:,np.isin(genes,markers)].mean()) if mask.any() else np.nan
        rows.append(r)
    return rows

def evaluate(scores,manifest):
    outputs=[]
    programs=list(PANELS)
    arms=sorted(manifest.treatment_arm.unique())
    for compartment in ['all_qc','cd8_like']:
        sub=scores[(scores.compartment==compartment)&scores.eligible]
        counts=sub.groupby('patient_id').timepoint.nunique()
        ids=sorted(counts[counts==2].index)
        if len(ids)<5:
            continue
        pre=sub[sub.timepoint=='pre'].set_index('patient_id').loc[ids,programs].to_numpy()
        post=sub[sub.timepoint=='post'].set_index('patient_id').loc[ids,programs].to_numpy()
        treatment=manifest.set_index('patient_id').loc[ids,'treatment_arm'].to_numpy()
        indicators=np.array([[int(t==a) for a in arms] for t in treatment])
        for i,patient in enumerate(ids):
            train=np.arange(len(ids))!=i
            assert patient not in np.array(ids)[train]
            assert set(treatment[train])==set(arms)
            delta=post[train]-pre[train]
            same=treatment[train]==treatment[i]
            assert same.any()
            scaler=StandardScaler().fit(pre[train])
            xtrain=np.column_stack([scaler.transform(pre[train]),indicators[train]])
            xtest=np.column_stack([scaler.transform(pre[[i]]),indicators[[i]]])
            ridge=Ridge(alpha=10).fit(xtrain,delta)
            predictions={'no_change':pre[i].copy(),
                         'post_mean':post[train].mean(axis=0),
                         'arm_post_mean':post[train][same].mean(axis=0),
                         'mean_delta':pre[i]+delta.mean(axis=0),
                         'arm_mean_delta':pre[i]+delta[same].mean(axis=0),
                         'ridge_delta':pre[i]+ridge.predict(xtest)[0]}
            scale=np.maximum(post[train].std(axis=0,ddof=1),1e-8)
            for model,prediction in predictions.items():
                assert np.isfinite(prediction).all()
                if model=='no_change':
                    np.testing.assert_array_equal(prediction,pre[i])
                for k,program in enumerate(programs):
                    outputs.append(dict(compartment=compartment,patient_id=patient,model=model,
                                        program=program,truth=post[i,k],prediction=prediction[k],
                                        absolute_error=abs(post[i,k]-prediction[k]),
                                        scaled_absolute_error=abs(post[i,k]-prediction[k])/scale[k],
                                        n_training_patients=int(train.sum())))
    predictions=pd.DataFrame(outputs)
    predictions.to_csv(OUT/'patient_heldout_predictions.tsv',sep='\t',index=False)
    summary=predictions.groupby(['compartment','model']).agg(
        n_patients=('patient_id','nunique'),mean_scaled_absolute_error=('scaled_absolute_error','mean'))
    summary.to_csv(OUT/'baseline_metrics.tsv',sep='\t')
    predictions.groupby(['compartment','model','program']).absolute_error.mean().to_csv(
        OUT/'per_program_mae.tsv',sep='\t')
    comparisons=[]
    rng=np.random.default_rng(20260909)
    for compartment,group in predictions.groupby('compartment'):
        errors=group.groupby(['patient_id','model']).scaled_absolute_error.mean().unstack()
        for baseline in ['no_change','mean_delta','arm_mean_delta','post_mean','arm_post_mean']:
            difference=(errors.ridge_delta-errors[baseline]).to_numpy()
            boot=rng.choice(difference,size=(10000,len(difference)),replace=True).mean(axis=1)
            comparisons.append(dict(compartment=compartment,comparison='ridge_delta minus '+baseline,
                                    mean_error_difference=float(difference.mean()),
                                    bootstrap_lower=float(np.quantile(boot,.025)),
                                    bootstrap_upper=float(np.quantile(boot,.975)),
                                    patients_improved=int((difference<0).sum()),n_patients=len(difference),
                                    uncertainty='conditional_on_LOOCV_fits_not_external_validation'))
    pd.DataFrame(comparisons).to_csv(OUT/'paired_error_comparisons.tsv',sep='\t',index=False)
    print(summary.to_string(),flush=True)
    print(pd.DataFrame(comparisons).to_string(index=False),flush=True)

def main():
    RAW.mkdir(parents=True,exist_ok=True)
    manifest=pd.read_csv(OUT/'patient_manifest.tsv',sep='\t')
    downloads=pd.read_csv(OUT/'download_manifest.tsv',sep='\t').to_dict('records')
    downloaded,failures=[],[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures={pool.submit(fetch,r):r for r in downloads}
        for future in as_completed(futures):
            row=futures[future]
            try:
                downloaded.append(future.result())
            except Exception as exc:
                failures.append(dict(url=row['url'],error=str(exc)))
            if (len(downloaded)+len(failures))%10==0:
                print(f'files complete={len(downloaded)} failed={len(failures)} total={len(downloads)}',flush=True)
    pd.DataFrame(downloaded).to_csv(OUT/'download_checksums.tsv',sep='\t',index=False)
    (OUT/'download_failures.json').write_text(json.dumps(failures,indent=2),encoding='utf-8')
    if failures:
        raise RuntimeError('Download incomplete; stop before cohort analysis. Rerun is resumable.')
    scores=[]
    for patient in manifest.loc[manifest.candidate_pair,'patient_id']:
        for timepoint in ['pre','post']:
            scores.extend(score(patient,timepoint))
        print(f'scored {patient}',flush=True)
    scores=pd.DataFrame(scores)
    scores.to_csv(OUT/'sample_qc_program_scores.tsv',sep='\t',index=False)
    evaluate(scores,manifest)

if __name__=='__main__':
    main()

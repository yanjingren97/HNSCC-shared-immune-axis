"""Leakage-controlled exploratory LOPO incremental MAE, full bootstrap refitting."""
from pathlib import Path
import json, argparse
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/signature_revision_20260915/incremental'
SRC=ROOT/'results/published_signature_external_audit'
CFG=json.loads((ROOT/'config/signature_composition_revision_20260915.json').read_text(encoding='utf-8'))
SIGS=pd.read_csv(SRC/'signature_catalog.tsv',sep='\t').signature.tolist()
TARGETS=CFG['incremental_targets']

def transform(train,test):
    r=rankdata(train,axis=0)
    rt=np.array([(np.sum(train[:,j]<test[j])+.5*np.sum(train[:,j]==test[j])+.5) for j in range(train.shape[1])])
    sd=r.std(axis=0,ddof=1); sd=np.where(sd>0,sd,1)
    return (r-r.mean(axis=0))/sd,(rt-r.mean(axis=0))/sd

def predict(tr,te,y):
    design=np.column_stack([np.ones(len(y)),tr])
    coef=np.linalg.lstsq(design,y,rcond=None)[0]
    return float(np.clip(np.r_[1,te]@coef,0,100))

def crossvalidate(x,y,indices=None):
    if indices is None: indices=np.arange(len(y))
    unique,weights=np.unique(indices,return_counts=True)
    if len(unique)<5: return None,None
    records=[]
    for i,w in zip(unique,weights):
        train=indices[indices!=i] # exclude all replicated copies of held-out donor
        z,t=transform(x[train],x[i]); yy=y[train]
        one=dict(patient_index=i,weight=w,observed=y[i],intercept=np.clip(yy.mean(),0,100))
        cyt=SIGS.index('Cytotoxic2')
        one['Cytotoxic2']=predict(z[:,[cyt]],t[[cyt]],yy)
        for name in TARGETS:
            j=SIGS.index(name); other=[k for k in range(len(SIGS)) if k!=j]
            _,_,vt=np.linalg.svd(z[:,other],full_matrices=False)
            loading=vt[0]
            if loading.sum()<0: loading=-loading
            axis=z[:,other]@loading; axist=t[other]@loading
            one[f'Cytotoxic2_plus_{name}']=predict(z[:,[cyt,j]],t[[cyt,j]],yy)
            one[f'PC1_without_{name}']=predict(axis[:,None],np.asarray([axist]),yy)
            one[f'PC1_plus_{name}']=predict(np.column_stack([axis,z[:,j]]),np.asarray([axist,t[j]]),yy)
        records.append(one)
    frame=pd.DataFrame(records)
    models=[c for c in frame if c not in ['patient_index','weight','observed']]
    mae={m:float(np.average(np.abs(frame[m]-frame.observed),weights=frame.weight)) for m in models}
    return mae,frame

def main(nboot):
    OUT.mkdir(parents=True,exist_ok=True); rng=np.random.default_rng(CFG['seed'])
    legacy=pd.read_csv(SRC/'multicohort_ssgsea_scores.tsv',sep='\t')
    main=pd.read_csv(ROOT/'results/signature_revision_20260915/composition/scores_T_clean.tsv',sep='\t').set_index('patient_id')
    summaries=[]; predictions=[]; bootrows=[]
    for cohort in ['GSE288199','GSE286827','GSE296954']:
        frame=legacy[legacy.cohort.eq(cohort)].set_index('sample_key')
        x=(main.loc[frame.index,SIGS] if cohort=='GSE288199' else frame[SIGS]).to_numpy(float)
        y=frame.continuous_outcome.to_numpy(float)
        assert np.isfinite(x).all() and np.isfinite(y).all() and frame.index.is_unique
        mae,pred=crossvalidate(x,y); pred['patient_id']=frame.index[pred.patient_index]; pred['cohort']=cohort
        predictions.append(pred)
        comps=[]
        for target in TARGETS:
            if target=='GZMK40' and cohort=='GSE296954': continue
            comps.extend([(target,'Cytotoxic2',f'Cytotoxic2_plus_{target}'),
                          (target,f'PC1_without_{target}',f'PC1_plus_{target}')])
        bs={c:[] for c in comps}; valid=0
        for b in range(nboot):
            indices=rng.integers(0,len(y),len(y)); bm,_=crossvalidate(x,y,indices)
            if bm is None: continue
            valid+=1
            for c in comps: bs[c].append(bm[c[1]]-bm[c[2]])
            if (b+1)%250==0: print(cohort,'bootstrap',b+1,flush=True)
        for target,base,aug in comps:
            vals=bs[(target,base,aug)]; lo,hi=np.quantile(vals,[.025,.975])
            summaries.append(dict(cohort=cohort,n=len(y),target=target,baseline=base,augmented=aug,
                baseline_MAE=mae[base],augmented_MAE=mae[aug],improvement_MAE=mae[base]-mae[aug],
                ci_low=lo,ci_high=hi,intercept_MAE=mae['intercept'],n_boot_valid=valid,n_boot_attempted=nboot,
                status='exploratory; source cohort excluded for GZMK40',
                score='clean T ssGSEA ES' if cohort=='GSE288199' else 'existing ssGSEA; train-rank transform removes common NES scale'))
            bootrows.extend(dict(cohort=cohort,target=target,baseline=base,bootstrap=b,improvement_MAE=v) for b,v in enumerate(vals))
        print(cohort,mae,flush=True)
        pd.DataFrame(summaries).to_csv(OUT/'incremental_MAE.tsv',sep='\t',index=False)
    pd.concat(predictions).to_csv(OUT/'heldout_predictions.tsv',sep='\t',index=False)
    pd.DataFrame(bootrows).to_csv(OUT/'bootstrap_improvements.tsv.gz',sep='\t',index=False)
    (OUT/'audit.json').write_text(json.dumps(dict(seed=CFG['seed'],bootstrap=nboot,
        resampling='patients with replacement; all copies of test patient excluded from train; full empirical-rank/PCA/regression refit',
        clipping='0 to 100 percent',endpoint='continuous pathological response; no mixed-endpoint pooling',
        interpretation='Exploratory percentile intervals for CV performance difference; no independent validation; no equivalence claim.'),indent=2),encoding='utf-8')

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--bootstrap',type=int,default=2000)
    args=parser.parse_args()
    with threadpool_limits(limits=1): main(args.bootstrap)

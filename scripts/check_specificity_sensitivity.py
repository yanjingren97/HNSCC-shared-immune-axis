from specificity_common import *
from scipy.stats import t as tdist
rng=np.random.default_rng(20260921)
x=load_expr('CPTAC');eng=FastSSGSEA(x);z=eng.score(SETS)
r=pd.read_csv(R/'results/independent_validation_20260915/CPTAC_patient_references.tsv',sep='\t',index_col=0).loc[x.columns]
prot=pd.read_csv(R/'data/external/signature_validation_20260915/HS_CPTAC_HNSCC_Proteomics_TMT_Gene_level_Tumor.cct',sep='\t',index_col=0)
panel=prot.loc[['PSMB9','TAP1','TAP2','GBP1','IFIT1'],x.columns].T.values
rows=[]
for j,target in enumerate(NAMES[:8]):
    other={k:[g for g in gs if g not in SETS[target]] for k,gs in SETS.items() if k not in [target,SHARED]}
    other={k:g for k,g in other.items() if len(eng.ids(g))>=2};zz=eng.score(other)
    def effect(ix):
        p=panel[ix];state=((p-p.mean(axis=0))/p.std(axis=0,ddof=1)).mean(axis=1)
        c=np.column_stack([zz[ix],r.CD3_IHC.values[ix],r.PTPRC_protein.values[ix]])
        ok=np.isfinite(c).all(axis=1);a=rankdata(z[ix,j][ok]);b=rankdata(state[ok]);c=rankdata(c[ok],axis=0)
        m=np.column_stack([np.ones(len(a)),c]);ar=a-m@np.linalg.lstsq(m,a,rcond=None)[0];br=b-m@np.linalg.lstsq(m,b,rcond=None)[0]
        value=np.corrcoef(ar,br)[0,1];df=len(a)-np.linalg.matrix_rank(m)-1
        return value,2*tdist.sf(abs(value)*np.sqrt(df/(1-value**2)),df),len(a)
    e,p,n=effect(np.arange(len(r)));boot=[effect(rng.integers(0,len(r),len(r)))[0] for _ in range(2000)]
    lo,hi=np.quantile(boot,[.025,.975]);rows.append(dict(target=ALIASES[target],n=n,partial_r=e,ci_low=lo,ci_high=hi,p=p,adjustment='All other target-gene-disjoint scores jointly plus CD3 and CD45; Shared84 excluded',status='Sensitivity specified after primary results; no target selection'))
df=pd.DataFrame(rows);order=np.argsort(df.p);q=np.empty(len(df));q[order]=np.minimum.accumulate((df.p.values[order]*len(df)/np.arange(1,len(df)+1))[::-1])[::-1].clip(0,1);df['q']=q
df.to_csv(O/'all_other_scores_sensitivity.tsv',sep='\t',index=False);print(df.to_string(index=False),flush=True)
# Verify the pre-existing constant baseline uses only the other patients.
held=pd.read_csv(R/'results/signature_revision_20260915/incremental/heldout_predictions.tsv',sep='\t');checks=[]
for cohort,d in held.groupby('cohort'):
    assert d.patient_id.is_unique
    expected=(d.observed.sum()-d.observed)/(len(d)-1);assert np.allclose(expected,d.intercept)
    for model in [c for c in d if c not in ['patient_index','weight','observed','patient_id','cohort']]:
        if d[model].notna().all():checks.append(dict(cohort=cohort,n=len(d),model=model,LOPO_MAE=abs(d.observed-d[model]).mean(),constant_LOPO_MAE=abs(d.observed-d.intercept).mean(),improvement_over_constant=abs(d.observed-d.intercept).mean()-abs(d.observed-d[model]).mean()))
pd.DataFrame(checks).to_csv(O/'clinical_constant_baseline_audit.tsv',sep='\t',index=False)

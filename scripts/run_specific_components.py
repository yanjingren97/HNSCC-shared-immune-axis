from specificity_common import *
from scipy.stats import t as tdist

rng=np.random.default_rng(CFG['seed']+1)
refs=pd.read_csv(R/'results/independent_validation_20260915/CPTAC_patient_references.tsv',sep='\t',index_col=0)
x=load_expr('CPTAC');eng=FastSSGSEA(x);scores=eng.score(SETS)
protein=pd.read_csv(R/'data/external/signature_validation_20260915/HS_CPTAC_HNSCC_Proteomics_TMT_Gene_level_Tumor.cct',sep='\t',index_col=0)
panel=protein.loc[['PSMB9','TAP1','TAP2','GBP1','IFIT1'],x.columns].T.values
refs=refs.loc[x.columns]

def partial(a,b,c):
    a=rankdata(a);b=rankdata(b);c=np.asarray(c)
    if c.ndim==1:c=c[:,None]
    m=np.column_stack([np.ones(len(a)),rankdata(c,axis=0)])
    ar=a-m@np.linalg.lstsq(m,a,rcond=None)[0];br=b-m@np.linalg.lstsq(m,b,rcond=None)[0]
    r=float(np.corrcoef(ar,br)[0,1]);df=len(a)-np.linalg.matrix_rank(m)-1
    p=float(2*tdist.sf(abs(r)*np.sqrt(df/max(1-r*r,1e-15)),df))
    return r,p

def bh(p):
    p=np.asarray(p);order=np.argsort(p);out=np.empty(len(p));out[order]=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(0,1);return out

rows=[];shares=[];coverage=[]
for j,target in enumerate(NAMES):
    for design in ['other_PC1','target_gene_disjoint_PC1']:
        other={k:g for k,g in SETS.items() if k!=target and k!=SHARED}
        if design.startswith('target_gene'):
            other={k:[g for g in gs if g not in SETS[target]] for k,gs in other.items()}
            other={k:gs for k,gs in other.items() if len(eng.ids(gs))>=2}
        z=eng.score(other)
        for k,gs in other.items():coverage.append(dict(target=ALIASES[target],design=design,reference_set=ALIASES[k],measured=len(eng.ids(gs))))
        _,axis=pc(z);yt=rankdata(scores[:,j]);fit=np.column_stack([np.ones(len(yt)),axis]);res=yt-fit@np.linalg.lstsq(fit,yt,rcond=None)[0]
        shares.append(dict(target=ALIASES[target],design=design,n=len(yt),reference_sets=len(other),in_sample_R2=1-(res@res)/np.sum((yt-yt.mean())**2),interpretation='Descriptive rank variance explained, not cross-validated prediction'))
        def evaluate(ix):
            _,ax=pc(z[ix]);p=panel[ix];state=((p-p.mean(axis=0))/p.std(axis=0,ddof=1)).mean(axis=1)
            rr=refs.iloc[ix];a=scores[ix,j]
            tests=[('CD3_IHC',rr.CD3_IHC.values,[]),('CD45_protein',rr.PTPRC_protein.values,[]),('IFN_protein',state,[]),('IFN_protein_adjusted_CD3_CD45',state,[rr.CD3_IHC.values,rr.PTPRC_protein.values])]
            out=[]
            for name,b,extra in tests:
                c=np.column_stack([ax]+extra);mask=np.isfinite(b)&np.isfinite(c).all(axis=1)
                effect,pval=partial(a[mask],b[mask],c[mask]);out.append((name,mask.sum(),effect,pval))
            return out
        obs=evaluate(np.arange(len(x.columns)));boots=[]
        for b in range(CFG['bootstrap_replicates']):boots.append([r[2] for r in evaluate(rng.integers(0,len(x.columns),len(x.columns)))])
        ci=np.quantile(boots,[.025,.975],axis=0)
        for t,(name,n,r,p) in enumerate(obs):rows.append(dict(target=ALIASES[target],design=design,reference=name,n=n,partial_r=r,ci_low=ci[0,t],ci_high=ci[1,t],p=p,primary=target!=SHARED))
        print(ALIASES[target],design,'done',flush=True)
res=pd.DataFrame(rows);res['q']=np.nan
# Separate families by reference and reference-construction design, across eight non-derived targets.
for _,ix in res[res.primary].groupby(['design','reference']).groups.items():res.loc[ix,'q']=bh(res.loc[ix,'p'])
res.to_csv(O/'specific_component_associations.tsv',sep='\t',index=False)
pd.DataFrame(shares).to_csv(O/'shared_rank_variance.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(O/'target_disjoint_reference_coverage.tsv',sep='\t',index=False)

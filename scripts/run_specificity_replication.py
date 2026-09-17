from specificity_common import *
from scipy.stats import spearmanr
from runpy import run_path
import hashlib
d=R/'data/external/specificity_20260917';rng=np.random.default_rng(20260919)
x=pd.read_csv(d/'TCGA_HNSC_HiSeqV2.gz',sep='\t',index_col=0)
assert x.index.is_unique
samples=sorted(c for c in x.columns if c[13:15]=='01');first={}
for s in samples:first.setdefault(s[:12],s)
m=pd.read_csv(d/'TCGA_leukocyte_methylation.tsv',sep='\t',header=None,names=['cancer','sample','leukocyte'])
m=m[(m.cancer=='HNSC')&(m['sample'].str[13:15]=='01')].copy();m['patient']=m['sample'].str[:12]
y=m.groupby('patient').leukocyte.mean();methylation_unique_n=len(y);ids=sorted(set(y.index)&set(first));x=x[[first[p] for p in ids]];x.columns=ids;y=y.loc[ids]
assert not x.isna().any().any() and x.columns.is_unique
eng=FastSSGSEA(x);cov=pd.DataFrame([dict(signature=ALIASES[k],measured=len(eng.ids(g)),total=len(g),coverage=len(eng.ids(g))/len(g)) for k,g in SETS.items()]);cov.to_csv(O/'TCGA_coverage.tsv',sep='\t',index=False)
assert cov.coverage.min()>=.8,cov
z=eng.score(SETS);pd.DataFrame(z,index=ids,columns=NAMES).to_csv(O/'TCGA_scores.tsv',sep='\t')
rows=[];plot=pd.DataFrame({'patient':ids,'leukocyte_methylation':y.values})
for label,jj in [('all9',list(range(9))),('without_Shared84',list(range(8)))]:
    v,axis=pc(z[:,jj]);r,p=spearmanr(axis,y);boot=[]
    for _ in range(2000):
        ix=rng.integers(0,len(ids),len(ids));_,a=pc(z[ix][:,jj]);boot.append(spearmanr(a,y.values[ix]).statistic)
    lo,hi=np.quantile(boot,[.025,.975]);rows.append(dict(design=label,n=len(ids),PC1_fraction=v,r=r,ci_low=lo,ci_high=hi,p=p));plot[label+'_PC1']=axis
res=pd.DataFrame(rows);order=np.argsort(res.p);q=np.empty(2);q[order]=np.minimum.accumulate((res.p.values[order]*2/np.arange(1,3))[::-1])[::-1].clip(0,1);res['q']=q
res.to_csv(O/'TCGA_methylation_replication.tsv',sep='\t',index=False);plot.to_csv(O/'TCGA_replication_plot_data.tsv',sep='\t',index=False)
audit={'RNA_primary_samples':len(samples),'RNA_unique_primary_patients':len(first),'methylation_primary_records':len(m),'methylation_unique_patients':methylation_unique_n,'paired_patients':len(ids),'duplicate_methylation_patients':m.patient[m.patient.duplicated()].tolist(),'all9_min_coverage':float(cov.coverage.min()),'RNA_source':'UCSC Xena TCGA.HNSC.sampleMap/HiSeqV2 2017-10-13','reference':'https://gdc.cancer.gov/about-data/publications/panimmune','scope':'Methylation-derived leukocyte abundance, not independent functional-state validation'}
(O/'TCGA_replication_audit.json').write_text(json.dumps(audit,indent=2));print(res.to_string(index=False),flush=True)

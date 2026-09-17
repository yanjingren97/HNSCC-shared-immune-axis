from specificity_common import *
from scipy.stats import t as tdist
O=R/'results/ifng6_stability_20260918'
d=R/'data/external/ifng6_stability_20260918'
assert (R/'config/independent_ifng6_20260918.json').exists()
f=d/'41467_2024_55203_MOESM7_ESM.xlsx'
gm=pd.read_csv(d/'gencode_v33_gene_symbols.tsv',sep='\t').set_index('ENSG')['gene'].to_dict()
x=pd.read_excel(f,sheet_name='2_Cohort_transcriptome',index_col=0)
x.index=[gm.get(str(g).split('.')[0],None) for g in x.index]
x=x.loc[x.index.notna()].groupby(level=0).mean();x.columns=x.columns.astype(str)
mp=pd.read_excel(f,sheet_name='3_DiffExp',header=2,usecols=[0,1,2])
pm={str(r.Protein).strip():gm.get(str(r.ENSG).split('.')[0]) for r in mp.itertuples()}
p=pd.read_excel(f,sheet_name='1_Cohort_proteome',index_col=0)
p.index=[pm.get(str(g).strip()) for g in p.index]
p=p.loc[p.index.notna()].groupby(level=0).mean();p.columns=p.columns.astype(str)
ids=x.columns.intersection(p.columns);x=x[ids];p=p[ids]
target=NAMES[3];other={k:[g for g in gs if g not in SETS[target]] for k,gs in SETS.items() if k not in [target,SHARED]}
cover=[dict(reference=ALIASES[k],total=len(gs),measured=sum(g in x.index for g in gs)) for k,gs in {target:SETS[target],**other}.items()]
pd.DataFrame(cover).to_csv(O/'independent_reference_coverage.tsv',sep='\t',index=False)
required=['PSMB9','TAP1','TAP2','GBP1','IFIT1','PTPRC','CD3E']
audit={'paired_n':len(ids),'RNA_gene_n':len(x),'protein_gene_n':len(p),'missing_proteins':[g for g in required if g not in p.index],'reference_coverage':cover}
if audit['missing_proteins'] or cover[0]['measured']!=6 or any(r['measured']/r['total']<.8 or r['measured']<2 for r in cover):
 audit['status']='Failed frozen feasibility gate; no association tested'
else:
 valid=p.loc[required].notna().all(axis=0)&x.notna().all(axis=0)
 audit['complete_n']=int(valid.sum())
 if valid.sum()<40:audit['status']='Failed complete-case gate; no association tested'
 else:
  x=x.loc[:,valid];p=p.loc[required,valid].T;engine=FastSSGSEA(x)
  score=engine.score({target:SETS[target]}).ravel();z=engine.score(other)
  cov=np.column_stack([z,p[['PTPRC','CD3E']].values]);pv=p.iloc[:,:5].values;n=len(x.columns)
  def calc(ix):
   pp=pv[ix];y=((pp-pp.mean(0))/pp.std(0,ddof=1)).mean(1)
   a=rankdata(score[ix]);b=rankdata(y);m=np.column_stack([np.ones(len(ix)),standardized(cov[ix])]);rk=np.linalg.matrix_rank(m);df=len(ix)-rk-1
   a-=m@np.linalg.lstsq(m,a,rcond=None)[0];b-=m@np.linalg.lstsq(m,b,rcond=None)[0];rr=float(np.corrcoef(a,b)[0,1])
   return dict(r=rr,df=int(df),rank=int(rk),condition=float(np.linalg.cond(m)),p=float(2*tdist.sf(abs(rr)*np.sqrt(df/(1-rr**2)),df)))
  result=calc(np.arange(n));rng=np.random.default_rng(20260922);bs=[calc(rng.integers(0,n,n))['r'] for _ in range(2000)];result['ci_low'],result['ci_high']=map(float,np.quantile(bs,[.025,.975]));audit.update(status='Frozen primary analysis completed',result=result)
  pd.DataFrame({'patient':x.columns,'IFNG6':score}).to_csv(O/'independent_patient_scores_local.tsv',sep='\t',index=False)
(O/'independent_cohort_audit.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit,indent=2))

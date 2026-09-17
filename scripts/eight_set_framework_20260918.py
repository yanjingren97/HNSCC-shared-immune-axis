from specificity_common import *
from scipy.stats import spearmanr
from collections import Counter
O=R/'results/ifng6_stability_20260918';rng=np.random.default_rng(20260923);sets={k:SETS[k] for k in NAMES[:8]}
rows=[];null=[];summ=[];coverage=[];axes={}
for cohort in CFG['core_cohorts']:
 x=load_expr(cohort);eng=FastSSGSEA(x);z=eng.score(sets);variance,axis=pc(z);axes[cohort]=z
 pd.DataFrame(z,index=x.columns,columns=list(sets)).to_csv(O/(cohort+'_eight_scores.tsv'),sep='\t')
 pd.DataFrame(np.corrcoef(standardized(z),rowvar=False),index=[ALIASES[k] for k in sets],columns=[ALIASES[k] for k in sets]).to_csv(O/(cohort+'_eight_correlations.tsv'),sep='\t')
 if cohort=='CPTAC':
  r=pd.read_csv(R/'results/independent_validation_20260915/CPTAC_patient_references.tsv',sep='\t',index_col=0).loc[x.columns];r['PC1_eight']=axis;r.to_csv(O/'CPTAC_eight_references.tsv',sep='\t')
  for ref in ['CD3_IHC','PTPRC_protein']:
   keep=r[ref].notna().values;effect,p=spearmanr(axis[keep],r.loc[keep,ref]);bs=[]
   for _ in range(2000):
    ix=rng.integers(0,len(x.columns),len(x.columns));_,a=pc(z[ix]);ok=keep[ix];bs.append(spearmanr(a[ok],r[ref].values[ix][ok]).statistic)
   lo,hi=np.quantile(bs,[.025,.975]);rows.append(dict(cohort=cohort,reference=ref,n=int(keep.sum()),PC1_variance=variance,r=effect,ci_low=lo,ci_high=hi,p=p))
 occ=Counter(g for gs in sets.values() for g in gs)
 reduced={k:[g for g in gs if occ[g]==1 and g in eng.lookup] for k,gs in sets.items()};reduced={k:gs for k,gs in reduced.items() if len(gs)>=2}
 for k,gs in reduced.items():coverage.append(dict(cohort=cohort,signature=ALIASES[k],full_n=len(eng.ids(sets[k])),exclusive_n=len(gs),exclusive_genes=';'.join(gs)))
 full=pc(eng.score({k:sets[k] for k in reduced}))[0];observed=pc(eng.score(reduced))[0];v=[]
 for b in range(1000):
  zz=np.column_stack([eng.score_ids(rng.choice(eng.ids(sets[k]),len(gs),replace=False)) for k,gs in reduced.items()]);val=pc(zz)[0];v.append(val);null.append(dict(cohort=cohort,replicate=b,PC1=val))
 summ.append(dict(cohort=cohort,sets=len(reduced),full_PC1=full,exclusive_PC1=observed,null_median=np.median(v),null_low=np.quantile(v,.025),null_high=np.quantile(v,.975),p=(1+sum(a<=observed for a in v))/1001))
print('Eight-set references',rows,flush=True)
pd.DataFrame(rows).to_csv(O/'eight_set_abundance_associations.tsv',sep='\t',index=False);pd.DataFrame(summ).to_csv(O/'eight_set_deletion_summary.tsv',sep='\t',index=False);pd.DataFrame(null).to_csv(O/'eight_set_deletion_replicates.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(O/'eight_set_exclusive_coverage.tsv',sep='\t',index=False)
a=axes['GSE288199_native22'];b=axes['GSE288199_fixed22'];vals=[]
for _ in range(2000):
 ix=rng.integers(0,len(a),len(a));vals.append((pc(a[ix])[0],pc(b[ix])[0]))
vals=np.array(vals);delta=vals[:,1]-vals[:,0];out={'n':len(a),'native_PC1':pc(a)[0],'fixed_PC1':pc(b)[0],'paired_change':pc(b)[0]-pc(a)[0],'change_low':np.quantile(delta,.025),'change_high':np.quantile(delta,.975),'native_low':np.quantile(vals[:,0],.025),'native_high':np.quantile(vals[:,0],.975),'fixed_low':np.quantile(vals[:,1],.025),'fixed_high':np.quantile(vals[:,1],.975)}
(O/'eight_set_composition.json').write_text(json.dumps(out,indent=2));print(out,flush=True)

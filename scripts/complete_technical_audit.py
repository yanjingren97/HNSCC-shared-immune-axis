from specificity_common import *
from scipy.stats import mannwhitneyu,spearmanr
# Selection audit uses the existing patient-level observed endpoint solely as an audit.
t=pd.read_csv(O/'single_cell_technical_covariates.tsv',sep='\t');h=pd.read_csv(R/'results/signature_revision_20260915/incremental/heldout_predictions.tsv',sep='\t');h=h[h.cohort=='GSE288199'];t=t.merge(h[['patient_id','observed']],on='patient_id',validate='one_to_one')
a=t.loc[t.included,'observed'];b=t.loc[~t.included,'observed'];s=pd.read_csv(O/'selection_22_vs_5.tsv',sep='\t');s=pd.concat([s,pd.DataFrame([dict(variable='pathological_response_percent',n_included=len(a),n_excluded=len(b),included_median=a.median(),excluded_median=b.median(),p=mannwhitneyu(a,b).pvalue)])],ignore_index=True)
def bh(p):
 p=np.asarray(p);i=np.argsort(p);q=np.empty(len(p));q[i]=np.minimum.accumulate((p[i]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(0,1);return q
s['q']=bh(s.p);s.to_csv(O/'selection_22_vs_5.tsv',sep='\t',index=False)
rows=[];sub=[]
for cohort in ['GSE288199_native22','GSE288199_fixed22']:
 x=load_expr(cohort);z=standardized(FastSSGSEA(x).score(SETS));c=t.set_index('patient_id').loc[x.columns];_,axis=pc(z)
 for name,cols in [('UMI',['UMI_per_cell']),('UMI_and_genes',['UMI_per_cell','genes_per_cell']),('UMI_genes_cells',['UMI_per_cell','genes_per_cell','n_cells'])]:
  m=np.column_stack([np.ones(len(c)),rankdata(c[cols].values,axis=0)]);res=z-m@np.linalg.lstsq(m,z,rcond=None)[0]
  rows.append(dict(cohort=cohort,adjustment=name,n=len(c),unadjusted_PC1=pc(z)[0],residual_rank_PC1=pc(res)[0],note='Descriptive sensitivity; correlated technical measures can also reflect biology'))
 f=pd.read_csv(O/'residual_T_subgroup_fractions.tsv',sep='\t').pivot(index='patient_id',columns='label',values='fraction').fillna(0).loc[x.columns]
 for label in f:
  if (f[label]>0).sum()>=10:
   r,p=spearmanr(axis,f[label]);sub.append(dict(cohort=cohort,T_subgroup=label,n=len(c),r=r,p=p,reference='RNA-derived annotation; composition diagnostic, not orthogonal biological validation'))
pd.DataFrame(rows).to_csv(O/'depth_adjusted_shared_structure.tsv',sep='\t',index=False)
s=pd.DataFrame(sub);s['q']=np.nan
for _,ii in s.groupby('cohort').groups.items():s.loc[ii,'q']=bh(s.loc[ii,'p'])
s.to_csv(O/'residual_T_composition_associations.tsv',sep='\t',index=False)
s=pd.read_csv(O/'null_summary.tsv',sep='\t');s['q_all_30_tests']=bh(s.p);s.to_csv(O/'null_summary.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).to_string(index=False))

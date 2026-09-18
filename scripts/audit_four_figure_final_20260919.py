"""Verify six/eight-set baselines and unify technical sensitivity to eight scores."""
from specificity_common import *
from scipy.stats import spearmanr
from collections import Counter
import hashlib
Q=R/'results/ifng6_stability_20260918';A=R/'results/four_figure_final_audit_20260919';A.mkdir(exist_ok=True)
sets={k:SETS[k] for k in NAMES[:8]};occ=Counter(g for gs in sets.values() for g in gs)
previous=pd.read_csv(Q/'eight_set_deletion_summary.tsv',sep='\t').set_index('cohort')
rows=[];hashes=[]
def record(p):hashes.append({'file':p.relative_to(R).as_posix(),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
for cohort in CFG['core_cohorts']:
 x=load_expr(cohort);eng=FastSSGSEA(x);z=eng.score(sets)
 reduced={k:[g for g in gs if occ[g]==1 and g in eng.lookup] for k,gs in sets.items()};reduced={k:gs for k,gs in reduced.items() if len(gs)>=2}
 full8=pc(z)[0];full6=pc(eng.score({k:sets[k] for k in reduced}))[0];exclusive=pc(eng.score(reduced))[0]
 assert np.isclose(full6,previous.loc[cohort,'full_PC1'],atol=1e-12,rtol=0)
 assert np.isclose(exclusive,previous.loc[cohort,'exclusive_PC1'],atol=1e-12,rtol=0)
 stored=Q/(cohort+'_eight_scores.tsv');record(stored);scores=pd.read_csv(stored,sep='\t',index_col=0)
 assert np.allclose(z,scores.loc[x.columns,list(sets)].values,atol=1e-10,rtol=0)
 rows.append(dict(cohort=cohort,n_patients=len(x.columns),eight_set_PC1=full8,six_set_baseline_PC1=full6,six_set_exclusive_PC1=exclusive,six_set_names=';'.join(ALIASES[k] for k in reduced),excluded_set_names=';'.join(ALIASES[k] for k in sets if k not in reduced),difference_six_minus_eight_pp=(full6-full8)*100,baseline_recomputed_match=True))
 print(cohort,full8,full6,exclusive,flush=True)
pd.DataFrame(rows).to_csv(A/'six_vs_eight_PC1_audit.tsv',sep='\t',index=False)
deletion=previous.reset_index();order=np.argsort(deletion.p.values);adj=np.empty(len(deletion));adj[order]=np.minimum.accumulate((deletion.p.values[order]*len(deletion)/np.arange(1,len(deletion)+1))[::-1])[::-1].clip(0,1);deletion['q_three_tests']=adj
deletion.to_csv(A/'eight_set_deletion_summary_audited.tsv',sep='\t',index=False)
tpath=O/'single_cell_technical_covariates.tsv';record(tpath);tech=pd.read_csv(tpath,sep='\t').set_index('patient_id');out=[]
for cohort in ['GSE288199_native22','GSE288199_fixed22']:
 scores=pd.read_csv(Q/(cohort+'_eight_scores.tsv'),sep='\t',index_col=0);_,axis=pc(scores.values);t=tech.loc[scores.index]
 for col in ['n_cells','UMI_per_cell','genes_per_cell','mixed_fraction']:
  r,p=spearmanr(axis,t[col]);out.append(dict(cohort=cohort,framework='eight_non_derived_sets',covariate=col,n=len(t),r=r,p=p))
out=pd.DataFrame(out);order=np.argsort(out.p.values);q=np.minimum.accumulate((out.p.values[order]*len(out)/np.arange(1,len(out)+1))[::-1])[::-1];adj=np.empty(len(out));adj[order]=np.minimum(q,1);out['q_eight_technical_tests']=adj
out.to_csv(A/'eight_score_technical_associations.tsv',sep='\t',index=False)
report={'purpose':'Verify six-set deletion baseline from expression and make technical sensitivity use the existing eight-score axis; no target selection or new clinical modelling','inputs':hashes,'six_eight_audit':'six_vs_eight_PC1_audit.tsv','technical_associations':'eight_score_technical_associations.tsv','technical_BH_family':'Eight descriptive tests: four covariates in each of two nested mixtures; not independent replication','software':{'numpy':np.__version__,'pandas':pd.__version__},'analyses':'Exploratory eight-score technical recalculation performed during final consistency audit on 19 September 2026'}
(A/'audit_provenance.json').write_text(json.dumps(report,indent=2),encoding='utf8');print(out.to_string(index=False),flush=True)

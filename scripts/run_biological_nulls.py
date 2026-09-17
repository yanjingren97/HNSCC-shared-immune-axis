from specificity_common import *
from scipy.spatial.distance import cdist
from collections import Counter
import time,hashlib

rng=np.random.default_rng(CFG['seed']);B=CFG['random_replicates']
gmt=R/'data/external/specificity_20260917/GO_Biological_Process_2025.gmt'
import re
pattern=re.compile(r'immune|immunoglobulin|lymphocyte|leukocyte|antigen|interferon|cytokine|natural killer|myeloid|macrophage|dendritic|neutrophil|complement|\b[TB] Cell',re.I)
terms={line.split('\t')[0]:line.rstrip().split('\t')[2:] for line in gmt.read_text().splitlines() if pattern.search(line.split('\t')[0])}
immune=set().union(*map(set,terms.values()))
(O/'immune_background_definition.json').write_text(json.dumps({'source':'https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2025','sha256':hashlib.sha256(gmt.read_bytes()).hexdigest(),'rule':pattern.pattern,'terms':list(terms),'gene_count':len(immune),'note':'Broad GO immune-term union, not enrichment analysis. No score/reference association used to define background.'},indent=2))
pd.Series(sorted(immune)).to_csv(O/'immune_background_genes.tsv',sep='\t',index=False,header=['gene'])
allrows=[];summ=[];checks=[];matching=[];feas=[]
for cohort in CFG['core_cohorts']:
    start=time.time();x=load_expr(cohort);engine=FastSSGSEA(x)
    ids={k:engine.ids(g) for k,g in SETS.items()};z=engine.score(SETS)
    if cohort=='CPTAC':old=pd.read_csv(R/'results/independent_validation_20260915/CPTAC_scores.tsv',sep='\t',index_col=0).loc[x.columns,NAMES].values
    else:
        comp='native_four' if 'native' in cohort else 'fixed_equal_RNA_four'
        frame=pd.read_csv(R/'results/signature_revision_20260915/composition/all_compartment_scores.tsv',sep='\t')
        old=frame[(frame.analysis_set=='matched_four_50')&(frame.compartment==comp)].set_index('patient_id').loc[x.columns,NAMES].values
    checks.append(dict(cohort=cohort,max_ES_error=float(abs(z-old).max()),rank_error=float(abs(rankdata(z,axis=0)-rankdata(old,axis=0)).max())))
    assert abs(z-old).max()<1e-6
    occ=Counter(i for a in ids.values() for i in a)
    exclusive={k:np.array([i for i in a if occ[i]==1],dtype=int) for k,a in ids.items()}
    elig=[k for k in ids if len(exclusive[k])>=2]
    feasible=[k for k in elig if len(ids[k])-len(exclusive[k])<=len(exclusive[k])]
    for k in NAMES:feas.append(dict(cohort=cohort,signature=ALIASES[k],measured=len(ids[k]),exclusive=len(exclusive[k]),delete_count=len(ids[k])-len(exclusive[k]),nonoverlap_delete_feasible=k in feasible))
    designs={'all9':list(range(9)),'without_Shared84':list(range(8))}
    full={key:pc(z[:,idx])[0] for key,idx in designs.items()}
    for keys,tag,nononly in [(elig,'random_delete_all',False),(feasible,'random_delete_exclusive_only',True)]:
        observed=pc(np.column_stack([engine.score_ids(exclusive[k]) for k in keys]))[0]
        fullpc=pc(np.column_stack([engine.score_ids(ids[k]) for k in keys]))[0]
        vals=[]
        for b in range(B):
            sampled=[]
            for k in keys:
                if nononly:
                    remove=set(rng.choice(exclusive[k],len(ids[k])-len(exclusive[k]),replace=False));keep=np.array([i for i in ids[k] if i not in remove])
                else:keep=rng.choice(ids[k],len(exclusive[k]),replace=False)
                sampled.append(engine.score_ids(keep))
            value=pc(np.column_stack(sampled))[0];vals.append(value)
            allrows.append(dict(cohort=cohort,control=tag,design=f'{len(keys)}_same_sets',replicate=b,pc1=value))
        summ.append(dict(cohort=cohort,control=tag,design=f'{len(keys)}_same_sets',observed_pc1=observed,full_pc1=fullpc,null_median=np.median(vals),null_lo=np.quantile(vals,.025),null_hi=np.quantile(vals,.975),p=(1+np.sum(np.array(vals)<=observed))/(B+1),tail='actual reduced PC1 below deletion null'))
    # Mean expression, expression SD and zero/nonzero detection are outcome-independent.
    # In CPTAC these are deposited-value proxies, not an assay detection limit.
    features=np.column_stack([rankdata(x.mean(axis=1))/len(x),rankdata(x.std(axis=1))/len(x),(x.values>0).mean(axis=1)])
    union=sorted(set().union(*map(set,ids.values())));variable=x.std(axis=1).values>0
    for bg in ['all_measurable','immune_GO']:
        pool=np.array([i for i,g in enumerate(x.index) if variable[i] and i not in union and (bg=='all_measurable' or g in immune)],dtype=int)
        assert len(pool)>len(union)
        distance=cdist(features[union],features[pool])
        near=np.argsort(distance,axis=1)[:,:100]
        candidates={u:pool[near[j]] for j,u in enumerate(union)}
        for architecture in ['independent_sets','preserved_overlap']:
            values={key:[] for key in designs};errors=[]
            for b in range(B):
                mapping={};used=set()
                if architecture=='preserved_overlap':
                    for u in rng.permutation(union):
                        available=[v for v in candidates[u] if v not in used]
                        v=int(rng.choice(available[:30]));mapping[u]=v;used.add(v)
                sampled=[]
                for k in NAMES:
                    if architecture=='independent_sets':
                        mapping={};used=set()
                        for u in rng.permutation(ids[k]):
                            v=int(rng.choice([v for v in candidates[u] if v not in used][:30]));mapping[u]=v;used.add(v)
                    target=np.array([mapping[u] for u in ids[k]])
                    errors.append(abs(features[target]-features[ids[k]]).mean(axis=0))
                    sampled.append(engine.score_ids(target))
                mat=np.column_stack(sampled)
                for key,idx in designs.items():
                    v=pc(mat[:,idx])[0];values[key].append(v)
                    allrows.append(dict(cohort=cohort,control=bg+'_'+architecture,design=key,replicate=b,pc1=v))
            matching.append(dict(cohort=cohort,background=bg,architecture=architecture,pool_n=len(pool),mean_abs_expression_rank_difference=float(np.mean(errors,axis=0)[0]),mean_abs_SD_rank_difference=float(np.mean(errors,axis=0)[1]),mean_abs_detection_difference=float(np.mean(errors,axis=0)[2])))
            for key,vals in values.items():
                summ.append(dict(cohort=cohort,control=bg+'_'+architecture,design=key,observed_pc1=full[key],full_pc1=full[key],null_median=np.median(vals),null_lo=np.quantile(vals,.025),null_hi=np.quantile(vals,.975),p=(1+np.sum(np.array(vals)>=full[key]))/(B+1),tail='observed PC1 above matched null'))
            print(cohort,bg,architecture,'done',round(time.time()-start,1),flush=True)
    for key,df in [('null_replicates',allrows),('null_summary',summ),('ssgsea_equivalence',checks),('matching_diagnostics',matching),('deletion_feasibility',feas)]:pd.DataFrame(df).to_csv(O/(key+'.tsv'),sep='\t',index=False)
    print('Completed',cohort,flush=True)

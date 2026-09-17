from specificity_common import *
from scipy.stats import t as tdist
O=R/'results/ifng6_stability_20260918';O.mkdir(exist_ok=True)
D=Path.home()/'Desktop/免疫表达轴/生物学主线改稿_20260918';D.mkdir(exist_ok=True)
plan={'date':'2026-09-18','status':'Follow-up diagnostics after 17 September results; no model selection','primary':'Original five-protein panel and joint adjustment for all seven other target-gene-disjoint non-derived scores plus CD3/CD45','diagnostics':['standardized rank design condition number and rank','all adjustment correlations','delete each of 83 complete patients and refit ranks and protein standardization','individual proteins and leave-one-protein-out composites with 2000 patient bootstrap resamples'],'BH':'five single-protein tests and five deletion tests are separate descriptive families','eight_set_framework':'Chosen after provenance audit and first extension; nine-set results retained as sensitivity'}
(O/'plan.json').write_text(json.dumps(plan,indent=2))
x=load_expr('CPTAC');eng=FastSSGSEA(x);target=NAMES[3];score=eng.score({target:SETS[target]}).ravel()
other={k:[g for g in gs if g not in SETS[target]] for k,gs in SETS.items() if k not in [target,SHARED]}
assert all(len(eng.ids(gs))>=2 for gs in other.values())
z=eng.score(other);ref=pd.read_csv(R/'results/independent_validation_20260915/CPTAC_patient_references.tsv',sep='\t',index_col=0).loc[x.columns]
proteins=['PSMB9','TAP1','TAP2','GBP1','IFIT1'];p=pd.read_csv(R/'data/external/signature_validation_20260915/HS_CPTAC_HNSCC_Proteomics_TMT_Gene_level_Tumor.cct',sep='\t',index_col=0).loc[proteins,x.columns].T.values
cov=np.column_stack([z,ref.CD3_IHC,ref.PTPRC_protein]);covnames=[ALIASES[k] for k in other]+['CD3_IHC','CD45_protein'];valid=np.isfinite(cov).all(axis=1);n=len(x.columns)
def calc(ix,keep=range(5)):
    pp=p[ix][:,list(keep)];yy=((pp-pp.mean(axis=0))/pp.std(axis=0,ddof=1)).mean(axis=1)
    ok=np.isfinite(cov[ix]).all(axis=1);a=rankdata(score[ix][ok]);b=rankdata(yy[ok]);c=standardized(cov[ix][ok]);m=np.column_stack([np.ones(len(a)),c]);rank=np.linalg.matrix_rank(m);df=len(a)-rank-1
    ar=a-m@np.linalg.lstsq(m,a,rcond=None)[0];br=b-m@np.linalg.lstsq(m,b,rcond=None)[0];r=np.corrcoef(ar,br)[0,1]
    return dict(r=r,p=2*tdist.sf(abs(r)*np.sqrt(df/(1-r*r)),df),n=len(a),rank=int(rank),df=int(df),condition=float(np.linalg.cond(m)))
base=calc(np.arange(n));(O/'strict_model_diagnostics.json').write_text(json.dumps(base,indent=2));assert base['n']==83 and base['rank']==10 and base['df']==72
c=standardized(cov[valid]);pd.DataFrame(c,index=x.columns[valid],columns=covnames).to_csv(O/'standardized_rank_adjustment_matrix.tsv',sep='\t');pd.DataFrame(np.corrcoef(c,rowvar=False),index=covnames,columns=covnames).to_csv(O/'adjustment_correlation_matrix.tsv',sep='\t')
coverage=[]
for j,(k,gs) in enumerate(other.items()):
    present=[g for g in gs if g in eng.lookup]
    coverage.append(dict(reference=ALIASES[k],original_catalog_n=len(SETS[k]),after_target_removal_n=len(gs),measured_after_removal_n=len(present),coverage_of_remaining=len(present)/len(gs),fraction_of_original=len(present)/len(SETS[k]),score_SD_83=z[valid,j].std(ddof=1),unique_score_values_83=len(np.unique(z[valid,j])),genes_with_nonzero_variance_83=int((x.loc[present,x.columns[valid]].std(axis=1)>0).sum()),retained_genes=';'.join(present)))
pd.DataFrame(coverage).to_csv(O/'strict_reference_coverage_variation.tsv',sep='\t',index=False)
loo=[]
for i in np.flatnonzero(valid):loo.append(dict(deleted_patient=x.columns[i],**calc(np.delete(np.arange(n),i))))
pd.DataFrame(loo).to_csv(O/'leave_one_patient_out.tsv',sep='\t',index=False)
rng=np.random.default_rng(20260918);draws=rng.integers(0,n,(2000,n));rows=[]
for typ,labels in [('single_protein',proteins),('leave_one_protein_out',proteins),('original_panel',['all_five'])]:
    for label in labels:
        keep=[proteins.index(label)] if typ=='single_protein' else ([j for j,g in enumerate(proteins) if g!=label] if typ=='leave_one_protein_out' else list(range(5)))
        result=calc(np.arange(n),keep);bs=np.array([calc(ix,keep)['r'] for ix in draws]);lo,hi=np.quantile(bs,[.025,.975]);rows.append(dict(design=typ,protein=label,**result,ci_low=lo,ci_high=hi));print(typ,label,result['r'],lo,hi,flush=True)
df=pd.DataFrame(rows);df['q']=np.nan
for typ,ix in df.groupby('design').groups.items():
    pv=df.loc[ix,'p'].values;order=np.argsort(pv);q=np.empty(len(pv));q[order]=np.minimum.accumulate((pv[order]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1].clip(0,1);df.loc[ix,'q']=q
df.to_csv(O/'protein_panel_sensitivity.tsv',sep='\t',index=False)
print('MODEL',base,'LOO RANGE',min(v['r'] for v in loo),max(v['r'] for v in loo),flush=True)

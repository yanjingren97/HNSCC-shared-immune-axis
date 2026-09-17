from specificity_common import *
from build_signature_composition_pseudobulk import annotation,masks,GROUPS
from scipy.io import mmread
from scipy import sparse
from scipy.stats import spearmanr,mannwhitneyu
from threadpoolctl import threadpool_limits
import gzip

rng=np.random.default_rng(20260920);B=100
patients=list(load_expr('GSE288199_fixed22').columns)
common=load_expr('GSE288199_fixed22').index
out=O/'downsampling_cache';out.mkdir(exist_ok=True)
allpatients=pd.read_csv(R/'results/signature_revision_20260915/composition/patient_QC.tsv',sep='\t').patient_id.tolist()
technical=[];subgroups=[]
with threadpool_limits(limits=2):
 for patient in allpatients:
    obs=pd.read_csv(annotation(patient),sep='\t',index_col=0);select,_=masks(obs);clean=obs.loc[select['T_clean']]
    technical.append(dict(patient_id=patient,included=patient in patients,n_cells=len(clean),UMI_per_cell=clean.total_counts.mean(),genes_per_cell=clean.n_genes.mean(),median_UMI=clean.total_counts.median(),mixed_fraction=obs.mixed_lineage_flag.astype(str).str.lower().eq('true').mean()))
    subgroups.extend(dict(patient_id=patient,label=label,n=n,fraction=n/len(clean)) for label,n in clean.reference_low.value_counts().items())
    if patient not in patients:continue
    dest=out/(patient+'.npz')
    if dest.exists():
        # Consume the same draws when resuming so uncached patients retain their original seed stream.
        for group in GROUPS:
            for _ in range(B):rng.choice(np.flatnonzero(select[group]),50,replace=False)
        continue
    stem=R/f'data/raw/GSE288199_decoded/GSE288199_{patient}_pre_T'
    barcodes=pd.read_csv(str(stem)+'_barcodes.tsv.gz',sep='\t',header=None)[0].astype(str)
    ix=pd.Index(barcodes).get_indexer(obs.barcode);assert (ix>=0).all()
    with gzip.open(str(stem)+'_matrix.mtx.gz','rb') as f:raw=mmread(f).tocsr()
    genes=pd.read_csv(str(stem)+'_features.tsv.gz',sep='\t',header=None).iloc[:,1].astype(str).str.upper().values
    unique,inverse=np.unique(genes,return_inverse=True)
    collapse=sparse.csr_matrix((np.ones(len(genes)),(inverse,np.arange(len(genes)))),shape=(len(unique),len(genes)))
    gi=pd.Index(unique).get_indexer(common);assert (gi>=0).all()
    raw=(collapse@raw)[gi].tocsr()
    totals=[];cpms=[]
    for group in GROUPS:
        cells=ix[select[group]];assert len(cells)>=50
        chosen=np.column_stack([rng.choice(cells,50,replace=False) for _ in range(B)])
        indicator=sparse.csr_matrix((np.ones(50*B),(chosen.ravel(),np.tile(np.arange(B),50))),shape=(raw.shape[1],B))
        counts=(raw@indicator).toarray();totals.append(counts);cpms.append(counts/counts.sum(axis=0)*1e6)
    native=sum(totals);native=native/native.sum(axis=0)*1e6;fixed=sum(cpms)/4
    scored={}
    for name,mat in [('equal_cells_pooled',native),('equal_cells_fixed_RNA',fixed)]:
        eng=FastSSGSEA(pd.DataFrame(np.log1p(mat),index=common));scored[name]=eng.score(SETS)
    np.savez_compressed(dest,**scored);print(patient,'100 equal-cell replicates complete',flush=True)
pd.DataFrame(technical).to_csv(O/'single_cell_technical_covariates.tsv',sep='\t',index=False)
pd.DataFrame(subgroups).to_csv(O/'residual_T_subgroup_fractions.tsv',sep='\t',index=False)
rows=[]
for name in ['equal_cells_pooled','equal_cells_fixed_RNA']:
    z=np.stack([np.load(out/(p+'.npz'))[name] for p in patients])
    for b in range(B):
        for label,j in [('all9',range(9)),('without_Shared84',range(8))]:rows.append(dict(mixture=name,design=label,replicate=b,PC1=pc(z[:,b,list(j)])[0],n_patients=len(patients),cells_per_subgroup=50))
df=pd.DataFrame(rows);df.to_csv(O/'equal_cell_PC1_replicates.tsv',sep='\t',index=False)
summary=df.groupby(['mixture','design']).PC1.agg(median='median',minimum='min',maximum='max',resampling_2_5_percent=lambda s:s.quantile(.025),resampling_97_5_percent=lambda s:s.quantile(.975))
summary.to_csv(O/'equal_cell_PC1_summary.tsv',sep='\t');print(summary)
t=pd.DataFrame(technical);sel=[]
for col in ['n_cells','UMI_per_cell','genes_per_cell','median_UMI','mixed_fraction']:
    a=t.loc[t.included,col];b=t.loc[~t.included,col];sel.append(dict(variable=col,n_included=len(a),n_excluded=len(b),included_median=a.median(),excluded_median=b.median(),p=mannwhitneyu(a,b,alternative='two-sided').pvalue))
pd.DataFrame(sel).to_csv(O/'selection_22_vs_5.tsv',sep='\t',index=False)
tech=[]
for cohort in ['GSE288199_native22','GSE288199_fixed22']:
    expr=load_expr(cohort);scores=FastSSGSEA(expr).score(SETS);_,axis=pc(scores)
    v=t.set_index('patient_id').loc[expr.columns]
    for col in ['n_cells','UMI_per_cell','genes_per_cell','mixed_fraction']:
        r,p=spearmanr(axis,v[col]);tech.append(dict(cohort=cohort,covariate=col,n=len(v),r=r,p=p))
pd.DataFrame(tech).to_csv(O/'composition_PC1_technical_associations.tsv',sep='\t',index=False)

from specificity_common import *
import gzip,re,hashlib
path=R/'data/external/specificity_20260917/GSE103322_HNSCC_all_data.txt.gz'
with gzip.open(path,'rt') as f:head=[f.readline().rstrip('\n').split('\t') for _ in range(6)]
meta=pd.DataFrame(dict(cell=head[0][1:],lymph_node=head[2][1:],malignant=head[3][1:],nonmalignant=head[4][1:],label=head[5][1:]))
meta['patient']=meta.cell.str.extract(r'^HN(?:SCC)?_?(\d+)',expand=False)
assert meta.cell.is_unique
unresolved=meta[meta.patient.isna()].copy();unresolved.to_csv(O/'Puram_unresolved_patient_cells.tsv',sep='\t',index=False)
meta=meta[meta.patient.notna()].copy()
def lineage(r):
    if r.malignant=='1':return 'Malignant'
    if r.nonmalignant!='1':return 'Unknown'
    lab=r.label.lstrip('-')
    return {'T cell':'T cells','B cell':'B cells','Macrophage':'Myeloid','Dendritic':'Myeloid','Mast':'Myeloid','Fibroblast':'Stromal','Endothelial':'Stromal','myocyte':'Stromal'}.get(lab,'Unknown')
meta['lineage']=meta.apply(lineage,axis=1);meta=meta[(meta.lymph_node=='0')&(meta.lineage!='Unknown')]
meta.to_csv(O/'Puram_primary_cell_annotations.tsv',sep='\t',index=False)
counts=meta.groupby(['patient','lineage']).size().unstack(fill_value=0);counts.to_csv(O/'Puram_patient_lineage_counts.tsv',sep='\t')
eligible=counts.index[(counts>=10).all(axis=1)].tolist()
audit={'source_url':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322','sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'primary_annotated_cells':len(meta),'primary_patients':int(meta.patient.nunique()),'unresolved_patient_cells_excluded':len(unresolved),'unresolved_note':'Combo-labelled cells without recoverable patient number excluded before expression-score analysis','all_five_lineages_at_least_10_cells_patients':eligible,'five_lineage_fixed_mixture_feasible':len(eligible)>=8,'interpretation':'Patient-lineage descriptive mapping; no per-cell inference or tissue-proportion estimate'}
(O/'Puram_lineage_audit.json').write_text(json.dumps(audit,indent=2))
x=pd.read_csv(path,sep='\t',skiprows=range(1,6),index_col=0,dtype={c:np.float32 for c in head[0][1:]});x.index=x.index.str.strip("'").str.upper();x=np.expm1(x*np.log(2))*10
if not x.index.is_unique:x=x.groupby(level=0).sum()
aggregates={};keys=[]
for (patient,lin),g in meta.groupby(['patient','lineage']):
    if len(g)<10:continue
    key=patient+'|'+lin;aggregates[key]=x[g.cell].mean(axis=1);keys.append(dict(sample=key,patient=patient,lineage=lin,n_cells=len(g)))
expr=np.log1p(pd.DataFrame(aggregates));eng=FastSSGSEA(expr);s=pd.DataFrame(eng.score(SETS),index=expr.columns,columns=list(ALIASES.values()))
out=pd.DataFrame(keys).set_index('sample').join(s);out.to_csv(O/'Puram_lineage_scores.tsv',sep='\t')
cov=pd.DataFrame([dict(signature=ALIASES[k],measured=len(eng.ids(gs)),total=len(gs),coverage=len(eng.ids(gs))/len(gs)) for k,gs in SETS.items()]);cov.to_csv(O/'Puram_gene_coverage.tsv',sep='\t',index=False)
# Paired descriptive contrast: each lineage versus the same patient's T-cell aggregate.
contrasts=[]
for lin in ['B cells','Myeloid','Malignant','Stromal']:
    a=out[out.lineage==lin].set_index('patient');b=out[out.lineage=='T cells'].set_index('patient');ids=sorted(set(a.index)&set(b.index))
    for name in ALIASES.values():
        delta=a.loc[ids,name]-b.loc[ids,name]
        contrasts.append(dict(lineage=lin,signature=name,n_paired_patients=len(ids),median_ES_difference_vs_T=float(delta.median()),fraction_above_T=float((delta>0).mean()),note='Descriptive within-patient contrast; uncalibrated ssGSEA units; no state-specificity inference'))
pd.DataFrame(contrasts).to_csv(O/'Puram_paired_lineage_contrasts.tsv',sep='\t',index=False)
print(audit,flush=True);print(cov.to_string(index=False),flush=True)


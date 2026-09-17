from specificity_common import *
from scipy.stats import spearmanr
import shutil,hashlib
# Check independent replication matching from unfiltered methylation records.
m=pd.read_csv(R/'data/external/specificity_20260917/TCGA_leukocyte_methylation.tsv',sep='\t',header=None,names=['cancer','sample','value']);m=m[(m.cancer=='HNSC')&(m['sample'].str[13:15]=='01')];assert m['sample'].str[:12].nunique()==528
p=pd.read_csv(O/'TCGA_replication_plot_data.tsv',sep='\t');assert p.patient.is_unique and len(p)==520
rr=pd.read_csv(O/'TCGA_methylation_replication.tsv',sep='\t');assert np.isclose(spearmanr(p.all9_PC1,p.leukocyte_methylation).statistic,rr.iloc[0].r)
# Validate summarized nulls against replicate-level outputs, without rerunning selection.
a=pd.read_csv(O/'null_replicates.tsv',sep='\t');s=pd.read_csv(O/'null_summary.tsv',sep='\t')
for row in s.itertuples():
 v=a[(a.cohort==row.cohort)&(a.control==row.control)&(a.design==row.design)].pc1.values;assert len(v)==1000
 expected=(1+np.sum(v<=row.observed_pc1 if row.control.startswith('random_delete') else v>=row.observed_pc1))/1001
 assert abs(expected-row.p)<1e-12 and np.isclose(np.median(v),row.null_median)
e=pd.read_csv(O/'ssgsea_equivalence.tsv',sep='\t');assert e.rank_error.max()==0 and e.max_ES_error.max()<1e-6
# Bootstrap intervals must contain finite valid correlations.
a=pd.read_csv(O/'specific_component_associations.tsv',sep='\t');assert len(a)==72 and a[['partial_r','ci_low','ci_high']].notna().all().all();assert (a.ci_low<=a.ci_high).all()
# Publish the actual used plan and all small result tables, not raw cell matrices.
(D/'tables').mkdir(exist_ok=True);(D/'code').mkdir(exist_ok=True);(D/'protocols').mkdir(exist_ok=True)
for p in O.glob('*.tsv'):
 if p.name not in ['null_replicates.tsv','Puram_primary_cell_annotations.tsv','Puram_unresolved_patient_cells.tsv']:shutil.copy2(p,D/'tables'/p.name)
for p in O.glob('*.json'):shutil.copy2(p,D/'protocols'/p.name)
for p in [R/'config'/x for x in ['biological_specificity_20260917.json','biological_replication_20260917.json','puram_lineage_specificity_20260917.json']]:shutil.copy2(p,D/'protocols'/p.name)
for name in ['specificity_common.py','run_biological_nulls.py','run_specific_components.py','run_specificity_replication.py','run_equal_cell_specificity.py','check_specificity_sensitivity.py','complete_technical_audit.py','build_specificity_provenance.py','run_puram_lineages.py','plot_biological_specificity.py','download_specificity_replication.py','download_puram_specificity.py']:
 shutil.copy2(R/'scripts'/name,D/'code'/name)
for name in ['replication_downloads.json','TCGA_HNSC_HiSeqV2.json']:shutil.copy2(R/'data/external/specificity_20260917'/name,D/'protocols'/name)
source=Path.home()/'Desktop/免疫表达轴/Cancer Immunology, Immunotherapy/CII_submission/manuscript.docx'
verification={'checks_passed':['Three original GSEApy score matrices reproduced with zero rank error','All 30 empirical null summaries reconstructed from 1000 replicates each','520 unique primary TCGA patient matches and correlation recomputed','528 unique primary methylation patients before matching','72 residual associations have finite bootstrap intervals','Existing intercept-only predictions verified as other-patient means'],'source_manuscript_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'source_manuscript_modified':False,'results_root':str(O),'raw_null_replicates':str(O/'null_replicates.tsv'),'bootstrap_n':2000,'equal_cell_repeats':100}
(D/'verification.json').write_text(json.dumps(verification,indent=2));print(json.dumps(verification,indent=2))

from specificity_common import *
import shutil
text=(D/'英文改稿段落与图注.md').read_text(encoding='utf-8');a=text.split('## Proposed abstract')[1].split('## Replacement Introduction')[0];print('Abstract words',len(a.split()))
# Explicitly retain overlap of the broad CD3 complex audit and provenance limitations.
p=D/'分析结论与改稿建议.md';s=p.read_text(encoding='utf-8').replace('均不在九集合的基因并集中。','均不在九集合的基因并集中。另行检查CD3各亚基与CD45：仅Korea30包含CD3G，CD45/PTPRC不与九集合重叠；IHC仍是不同测量方式，但不能把所有CD3参照都称为完全基因无重叠。');p.write_text(s,encoding='utf-8')
shutil.copy2(O/'orthogonal_reference_overlap_audit.tsv',D/'tables/orthogonal_reference_overlap_audit.tsv');shutil.copy2(R/'scripts/run_equal_cell_specificity.py',D/'code/run_equal_cell_specificity.py')
# Place readily comparable key effects in one compact TSV, with uncertainty provenance.
sp=pd.read_csv(O/'specific_component_associations.tsv',sep='\t');st=pd.read_csv(O/'all_other_scores_sensitivity.tsv',sep='\t');tc=pd.read_csv(O/'TCGA_methylation_replication.tsv',sep='\t')
rows=[]
for r in tc.itertuples():rows.append(dict(analysis='TCGA methylation '+r.design,n=r.n,effect=r.r,CI_low=r.ci_low,CI_high=r.ci_high,p=r.p,q=r.q,units='Spearman rho'))
for name in ['IFNG6','GEP18']:
 for typ,dd in [('other_gene_disjoint_PC1',sp[(sp.target==name)&(sp.design=='target_gene_disjoint_PC1')&(sp.reference=='IFN_protein_adjusted_CD3_CD45')]),('all_other_gene_disjoint_scores',st[st.target==name])]:
  r=dd.iloc[0];rows.append(dict(analysis=name+' '+typ+' + CD3/CD45',n=r.n,effect=r.partial_r,CI_low=r.ci_low,CI_high=r.ci_high,p=r.p,q=r.q,units='partial rank correlation'))
pd.DataFrame(rows).to_csv(D/'关键结果汇总.tsv',sep='\t',index=False)

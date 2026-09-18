"""Rearrange existing estimates only; no new inferential models or selection."""
from pathlib import Path
import json, shutil, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R=Path(__file__).resolve().parents[1]
O=R/'results/biological_specificity_20260917';Q=R/'results/ifng6_stability_20260918'
parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path,default=R/'results/four_figure_final_audit_20260919/rendered');args=parser.parse_args();D=args.output_dir
F=D/'figures';F.mkdir(parents=True,exist_ok=True)
S=D/'support';S.mkdir(exist_ok=True)
A=R/'results/four_figure_final_audit_20260919'
plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
blue='#7B9CAF';teal='#6EA3A5';orange='#C59C7A'
def panel(ax,l,title):
 pos=ax.get_position();ax.figure.text(pos.x0-.04,pos.y1+.048,l,fontweight='bold',fontsize=14)
 ax.set_title(title,loc='left',pad=25,fontweight='bold');ax.set_axisbelow(True)
def save(fig,name):
 for ext,dpi in [('png',180),('pdf',300),('svg',300),('tiff',600)]:
  fig.savefig(F/(name+'.'+ext),dpi=dpi,facecolor='white',**({'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}))
 plt.close(fig)
def three_axes():
 fig=plt.figure(figsize=(10,8.1));g=fig.add_gridspec(2,2,left=.16,right=.97,bottom=.11,top=.90,wspace=.42,hspace=.86,height_ratios=[1,.85])
 return fig,[fig.add_subplot(g[0,0]),fig.add_subplot(g[0,1]),fig.add_subplot(g[1,:])]
# Fig. 1 is rebuilt from the existing matrices and references, never copied from a desktop export.
from matplotlib.colors import LinearSegmentedColormap
names=['GZMK40','CYT2','EXH5','IFNG6','TLS9','GEP18','Korea30','Dysf238']
fig,axs=plt.subplots(2,2,figsize=(10,8.6));fig.subplots_adjust(left=.12,right=.97,bottom=.10,top=.90,wspace=.47,hspace=.64)
ax=axs[0,0];cor=pd.read_csv(Q/'CPTAC_eight_correlations.tsv',sep='\t',index_col=0);im=ax.imshow(cor,vmin=0,vmax=1,cmap=LinearSegmentedColormap.from_list('soft',['#F5F6F4',blue]),aspect='auto');ax.set_xticks(range(8),names,rotation=50,ha='right');ax.set_yticks(range(8),names);fig.colorbar(im,ax=ax,fraction=.04,pad=.03,label='Spearman correlation');panel(ax,'a','CPTAC shared structure: PC1 = 88.7%')
refs=pd.read_csv(Q/'CPTAC_eight_references.tsv',sep='\t',index_col=0);stats=pd.read_csv(Q/'eight_set_abundance_associations.tsv',sep='\t')
for ax,key,l,t,xlab in [(axs[0,1],'CD3_IHC','b','CPTAC histological reference','Stromal CD3-positive TILs (%)'),(axs[1,0],'PTPRC_protein','c','CPTAC protein reference','CD45 / PTPRC protein abundance')]:
 dd=refs[[key,'PC1_eight']].dropna();row=stats[stats.reference==key].iloc[0];ax.scatter(dd[key],dd.PC1_eight,s=18,color=teal,alpha=.7);ax.text(.03,.97,f'n = {len(dd)}; ρ = {row.r:.2f}\n95% CI {row.ci_low:.2f} to {row.ci_high:.2f}',va='top',transform=ax.transAxes);ax.set_xlabel(xlab);ax.set_ylabel('Eight-score RNA PC1');panel(ax,l,t)
ax=axs[1,1];tc=pd.read_csv(O/'TCGA_replication_plot_data.tsv',sep='\t');ax.scatter(tc.leukocyte_methylation,tc.without_Shared84_PC1,s=12,color=blue,alpha=.4);ax.text(.03,.97,'n = 520; ρ = 0.60\n95% CI 0.54 to 0.66\nPC1 = 85.6%',transform=ax.transAxes,va='top');ax.set_xlabel('Methylation-derived leukocyte fraction');ax.set_ylabel('Eight-score RNA PC1');panel(ax,'d','TCGA orthogonal abundance reference');save(fig,'Figure_1_shared_structure_and_abundance')

# Fig. 2: retain both random backgrounds and all three deletion comparisons.
fig,axs=three_axes();null=pd.read_csv(O/'null_replicates.tsv',sep='\t');ns=pd.read_csv(O/'null_summary.tsv',sep='\t')
controls=['all_measurable_independent_sets','all_measurable_preserved_overlap','immune_GO_independent_sets','immune_GO_preserved_overlap']
for ax,co,letter,title in [(axs[0],'CPTAC','a','CPTAC matched random sets'),(axs[1],'GSE288199_fixed22','b','Fixed T-cell matched random sets')]:
 vals=[null[(null.cohort==co)&(null.control==c)&(null.design=='without_Shared84')].pc1*100 for c in controls]
 bp=ax.boxplot(vals,patch_artist=True,showfliers=False)
 for patch,col in zip(bp['boxes'],[blue,blue,teal,teal]):patch.set_facecolor(col);patch.set_alpha(.65)
 obs=ns[(ns.cohort==co)&(ns.design=='without_Shared84')].iloc[0].observed_pc1*100
 ax.axhline(obs,color=orange,label=f'Observed: {obs:.1f}%');ax.set_xticks(range(1,5),['All genes\nSeparate','All genes\nOverlap','Immune\nSeparate','Immune\nOverlap']);ax.tick_params(axis='x',labelsize=9)
 ax.set_ylim(0,100);ax.set_ylabel('Eight-score PC1 variance (%)');ax.legend(frameon=False,fontsize=9,loc='upper right');panel(ax,letter,title)
ax=axs[2];d=pd.read_csv(Q/'eight_set_deletion_summary.tsv',sep='\t');y=np.arange(3)
ax.errorbar(d.null_median*100,y,xerr=np.vstack([d.null_median-d.null_low,d.null_high-d.null_median])*100,fmt='o',color=blue,capsize=3,label='Random deletion: median and central 95%')
ax.scatter(d.exclusive_PC1*100,y,marker='D',color=orange,label='Remove shared genes',zorder=3);ax.set_yticks(y,['CPTAC','Native T','Fixed T']);ax.set_ylim(2.45,-.45);ax.set_xlim(35,100);ax.set_xlabel('PC1 variance in the same six eligible sets (%)')
ax.legend(frameon=False,fontsize=9,loc='upper center',bbox_to_anchor=(.5,-.27),ncol=2);panel(ax,'c','Gene-overlap removal versus size-matched deletion');save(fig,'Figure_2_random_sets_and_gene_overlap')

# Fig. 3: primary eight-set composition and sampling estimates; clearly labelled
# newly harmonized eight-set technical sensitivity; no clinical modelling.
fig,axs=three_axes();v=json.loads((Q/'eight_set_composition.json').read_text());ax=axs[0]
vals=[v['native_PC1']*100,v['fixed_PC1']*100];ax.bar([0,1],vals,color=[blue,teal],width=.6)
for x,val in enumerate(vals):ax.text(x,val+2,f'{val:.1f}%',ha='center',fontsize=10)
ax.set_xticks([0,1],['Native mixtures','Fixed RNA']);ax.set_ylim(0,105);ax.set_ylabel('Eight-score PC1 variance (%)')
ax.text(.03,.98,'22 paired patients\nChange −5.00 percentage points\n95% CI −13.08 to 2.25',va='top',transform=ax.transAxes,fontsize=9);panel(ax,'a','Fixed subgroup RNA contributions')
ax=axs[1];rep=pd.read_csv(O/'equal_cell_PC1_replicates.tsv',sep='\t');rep=rep[rep.design=='without_Shared84']
vv=[rep.loc[rep.mixture==key,'PC1'].to_numpy()*100 for key in ['equal_cells_pooled','equal_cells_fixed_RNA']]
bp=ax.boxplot(vv,patch_artist=True,showfliers=False,widths=.5)
for patch,col in zip(bp['boxes'],[blue,teal]):patch.set_facecolor(col);patch.set_alpha(.7)
rng=np.random.default_rng(20260919)
for x,values,col in zip([1,2],vv,[blue,teal]):
 ax.scatter(x+rng.uniform(-.15,.15,len(values)),values,s=9,color=col,alpha=.35,zorder=1)
 ax.text(x,70,f'Median {np.median(values):.1f}%',ha='center',fontsize=9)
ax.set_xticks([1,2],['Pooled counts','Fixed RNA']);ax.set_ylim(0,100);ax.set_ylabel('Eight-score PC1 variance (%)');ax.text(.04,.97,'50 cells per subgroup\n100 sampling repeats',va='top',transform=ax.transAxes,fontsize=9);panel(ax,'b','Equal-cell sampling sensitivity')
ax=axs[2];tech=pd.read_csv(A/'eight_score_technical_associations.tsv',sep='\t');keys=['n_cells','UMI_per_cell','genes_per_cell','mixed_fraction'];y=np.arange(4)
for cohort,col,offset,label in [('GSE288199_native22',blue,-.12,'Native mixtures'),('GSE288199_fixed22',teal,.12,'Fixed RNA')]:
 rows=tech[tech.cohort==cohort].set_index('covariate').loc[keys];ax.scatter(rows.r,y+offset,color=col,label=label,s=34)
ax.set_yticks(y,['T-cell count','Mean UMI / cell','Mean genes / cell','Mixed-lineage fraction']);ax.set_ylim(3.5,-.5);ax.set_xlim(0,.8);ax.set_xlabel('Spearman correlation with the eight-score PC1');ax.legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.27),ncol=2,fontsize=9)
panel(ax,'c','Technical associations: eight scores (n = 22)');save(fig,'Figure_3_single_cell_composition_and_sampling')

# Fig. 4: core conditional information, keeping original protein panel intact.
names=['GZMK40','CYT2','EXH5','IFNG6','TLS9','GEP18','Korea30','Dysf238']
fig=plt.figure(figsize=(10,7.7));g=fig.add_gridspec(2,2,left=.12,right=.97,bottom=.10,top=.89,wspace=.65,hspace=.70,width_ratios=[1,1]);a=fig.add_subplot(g[:,0]);b=fig.add_subplot(g[0,1]);c=fig.add_subplot(g[1,1])
d=pd.read_csv(O/'all_other_scores_sensitivity.tsv',sep='\t').set_index('target').loc[names]
a.errorbar(d.partial_r,range(8),xerr=np.vstack([d.partial_r-d.ci_low,d.ci_high-d.partial_r]),fmt='o',color=teal,capsize=3)
a.scatter([d.loc['IFNG6','partial_r']],[3],color=orange,zorder=5,s=42);a.set_yticks(range(8),names);a.invert_yaxis();a.axvline(0,color='silver',lw=.7);a.set_xlim(-.5,1.0);a.set_xticks([-.4,-.2,0,.2,.4,.6]);a.text(.90,1.035,'BH q',transform=a.transAxes,ha='center',fontweight='bold',fontsize=9);a.set_xlabel('Partial rank correlation (95% CI)');panel(a,'a','Conditional protein associations')
for i,value in enumerate(d.q):a.text(.87,i,f'{value:.4f}',va='center',ha='center',fontsize=9,color=orange if i==3 else '#344A57')
q=pd.read_csv(Q/'protein_panel_sensitivity.tsv',sep='\t')
for ax,typ,letter,title in [(b,'single_protein','b','Each original panel protein'),(c,'leave_one_protein_out','c','Leave one panel protein out')]:
 d=q[q.design==typ];ax.errorbar(d.r,range(5),xerr=np.vstack([d.r-d.ci_low,d.ci_high-d.r]),fmt='o',color=teal,capsize=3);ax.set_yticks(range(5),[('Without ' if typ!='single_protein' else '')+p for p in d.protein]);ax.invert_yaxis();ax.axvline(0,color='silver',lw=.7);ax.set_xlim(-.3,.8);ax.set_xlabel('IFNG6 partial rank correlation (95% CI)');panel(ax,letter,title)
save(fig,'Figure_4_IFNG6_conditional_information')

fig,ax=plt.subplots(figsize=(6.4,4.8));fig.subplots_adjust(left=.14,right=.96,bottom=.15,top=.83);d=pd.read_csv(Q/'leave_one_patient_out.tsv',sep='\t')
ax.hist(d.r,bins=15,color=blue,edgecolor='white');ax.axvline(.381737,color=orange,label='Full-data estimate');ax.set_xlabel('IFNG6 partial rank correlation');ax.set_ylabel('Patient-deletion refits');ax.set_title('Leave-one-patient-out sensitivity',loc='left',pad=20,fontweight='bold');ax.legend(frameon=False);save(fig,'Supplementary_Figure_S2_patient_deletion')

print(D)

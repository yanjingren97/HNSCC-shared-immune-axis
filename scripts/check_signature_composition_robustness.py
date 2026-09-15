"""Gene-disjoint and technical-depth controls for shared signature structure."""
from pathlib import Path
import sys,json
import numpy as np
import pandas as pd
from scipy.stats import rankdata,spearmanr
from analyze_signature_composition import structure,SETS,NAMES,score,OUT,CFG
import analyze_signature_composition as core

def main():
    rng=np.random.default_rng(CFG['seed']+1)
    unique={k:sorted(set(v)-set().union(*(set(v2) for k2,v2 in SETS.items() if k2!=k))) for k,v in SETS.items()}
    eligible={k:v for k,v in unique.items() if len(v)>=2}
    pd.DataFrame([dict(signature=k,original_genes=len(SETS[k]),exclusive_genes=len(v),evaluable=len(v)>=2,genes=';'.join(v)) for k,v in unique.items()]).to_csv(OUT/'exclusive_gene_catalog.tsv',sep='\t',index=False)
    counts=pd.read_csv(OUT/'compartment_counts.tsv',sep='\t')
    allscores=pd.read_csv(OUT/'all_compartment_scores.tsv',sep='\t')
    rows=[]; depths=[]; loadings=[]
    for tag in ['matched_four_50','matched_four_100']:
        for comp in ['native_four','fixed_equal_RNA_four']:
            expr=pd.read_csv(OUT/f'{tag}_{comp}_logCPM.tsv.gz',sep='\t',index_col=0)
            ids=expr.columns.tolist()
            core.GENES=expr.index.to_numpy(); core.CURRENT_IDS=ids; core.SETS=eligible; core.NAMES=list(eligible)
            ds=score(np.expm1(expr.to_numpy()))
            full=allscores[(allscores.analysis_set==tag)&(allscores.compartment==comp)].set_index('patient_id').loc[ids,list(eligible)]
            for kind,frame in [('full_same_six',full),('exclusive_six',ds)]:
                m,pc,loading,_=structure(frame.to_numpy()); bs=[]
                for _ in range(2000):
                    idx=rng.integers(0,len(frame),len(frame)); bs.append(structure(frame.to_numpy()[idx])[0]['pc1_variance'])
                lo,hi=np.quantile(bs,[.025,.975])
                rows.append(dict(analysis_set=tag,compartment=comp,gene_design=kind,n=len(frame),n_signatures=len(eligible),pc1_ci_low=lo,pc1_ci_high=hi,**m))
                frame.rename_axis('patient_id').to_csv(OUT/f'{tag}_{comp}_{kind}_scores.tsv',sep='\t')
                loadings.extend(dict(analysis_set=tag,compartment=comp,gene_design=kind,signature=k,loading=v) for k,v in zip(frame.columns,loading))
            print(tag,comp,'full6',structure(full.to_numpy())[0],'exclusive6',structure(ds.to_numpy())[0],flush=True)
    for comp in ['T_clean','CD8_clean']+list(CFG['primary_groups']):
        frame=allscores[(allscores.analysis_set=='all_eligible_50')&(allscores.compartment==comp)].set_index('patient_id')[NAMES]
        meta=counts[counts.compartment==comp].set_index('patient_id').loc[frame.index]
        cells=np.log1p(meta.n_cells.to_numpy()); umi_per_cell=np.log1p(meta.total_UMI/meta.n_cells)
        cov=rankdata(np.column_stack([cells,umi_per_cell]),axis=0)
        design=np.column_stack([np.ones(len(frame)),cov])
        z=rankdata(frame.to_numpy(),axis=0); residual=z-design@np.linalg.lstsq(design,z,rcond=None)[0]
        before,pc,_,_=structure(frame.to_numpy()); after,*_=structure(residual)
        depths.append(dict(compartment=comp,n=len(frame),pc1_before=before['pc1_variance'],pc1_after_depth_residualization=after['pc1_variance'],
            rho_PC1_cell_count=spearmanr(pc,cells).statistic,rho_PC1_UMI_per_cell=spearmanr(pc,umi_per_cell).statistic))
    pd.DataFrame(rows).to_csv(OUT/'exclusive_gene_robustness.tsv',sep='\t',index=False)
    pd.DataFrame(depths).to_csv(OUT/'technical_depth_robustness.tsv',sep='\t',index=False)
    pd.DataFrame(loadings).to_csv(OUT/'exclusive_gene_loadings.tsv',sep='\t',index=False)
    (OUT/'robustness_audit.json').write_text(json.dumps(dict(exclusive_signature_count=len(eligible),
        excluded_because_fewer_than_two_unique_genes=[k for k in unique if k not in eligible],
        status='exploratory structural sensitivity; not validation of original algorithms or signature specificity',
        depth_adjustment='rank residuals on log cell count and log UMI per cell; descriptive, no causal interpretation'),indent=2),encoding='utf-8')

if __name__=='__main__': main()

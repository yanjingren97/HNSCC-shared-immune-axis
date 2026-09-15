"""Post hoc extension of the already-defined six disjoint gene sets to CPTAC."""
import json,shutil
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from analyze_independent_validation import ROOT,D,O,SETS,NAMES,gp,structure
def main():
    out=O/'gene_disjoint';out.mkdir(exist_ok=True)
    plan={'status':'exploratory sensitivity added after the original CPTAC associations','design':'Use identical six eligible globally disjoint sets as prior composition analysis; compare same six full sets; no response outcomes; paired bootstrap of PC1 variance difference and IHC correlation difference; 2000 resamples','seed':20260917}
    (out/'plan.json').write_text(json.dumps(plan,indent=2))
    cat=pd.read_csv(ROOT/'results/signature_revision_20260915/composition/exclusive_gene_catalog.tsv',sep='\t')
    cat=cat[cat.evaluable];names=cat.signature.tolist();sets={r.signature:r.genes.split(';') for r in cat.itertuples()}
    assert len(sets)==6
    assert sum(map(len,sets.values()))==len(set().union(*map(set,sets.values())))
    full=pd.read_csv(O/'CPTAC_scores.tsv',sep='\t',index_col=0)[names]
    rna=pd.read_csv(D/'HS_CPTAC_HNSCC_RNAseq_RSEM_UQ_log2_Tumor.cct',sep='\t',index_col=0).loc[:,full.index]
    coverage=pd.DataFrame([dict(signature=k,unique_genes=len(v),measured_genes=len(set(v)&set(rna.index)),genes=';'.join(v)) for k,v in sets.items()])
    assert coverage.measured_genes.min()>=2
    coverage.to_csv(out/'gene_catalog_and_coverage.tsv',sep='\t',index=False)
    result=gp.ssgsea(data=rna.copy().reset_index(names='Gene'),gene_sets=sets,outdir=None,no_plot=True,sample_norm_method='rank',min_size=2,max_size=500,permutation_num=0,threads=1,seed=plan['seed'],verbose=False)
    exclusive=result.res2d.pivot(index='Name',columns='Term',values='ES').astype(float).loc[full.index,names]
    exclusive.to_csv(out/'exclusive_six_scores.tsv',sep='\t')
    ref=pd.read_csv(O/'CPTAC_patient_references.tsv',sep='\t',index_col=0).loc[full.index,'CD3_IHC']
    frames={'full_same_six':full,'exclusive_six':exclusive};values={};loads={};rows=[]
    for name,frame in frames.items():
        m,pc,loading,_=structure(frame.values);values[name]=[m['pc1_variance'],spearmanr(pc[ref.notna()],ref.dropna()).statistic];loads[name]=loading
        rows.append(dict(design=name,n=len(frame),n_IHC=ref.notna().sum(),PC1_variance=values[name][0],IHC_rho=values[name][1]))
    rng=np.random.default_rng(plan['seed']);boot={k:[] for k in frames}
    for _ in range(2000):
        ix=rng.integers(0,len(full),len(full));br=ref.iloc[ix];valid=br.notna().values
        for name,frame in frames.items():
            m,pc,_,_=structure(frame.values[ix],ref=loads[name]);boot[name].append([m['pc1_variance'],spearmanr(pc[valid],br.iloc[np.flatnonzero(valid)]).statistic])
    for row in rows:
        lo,hi=np.quantile(boot[row['design']],[.025,.975],axis=0)
        row.update(PC1_ci_low=lo[0],PC1_ci_high=hi[0],IHC_ci_low=lo[1],IHC_ci_high=hi[1])
    diff=np.asarray(boot['exclusive_six'])-np.asarray(boot['full_same_six']);lo,hi=np.quantile(diff,[.025,.975],axis=0);point=np.array(values['exclusive_six'])-np.array(values['full_same_six'])
    delta=pd.DataFrame({'metric':['PC1_variance','IHC_rho'],'exclusive_minus_full':point,'ci_low':lo,'ci_high':hi})
    pd.DataFrame(rows).to_csv(out/'structure_and_IHC.tsv',sep='\t',index=False);delta.to_csv(out/'paired_differences.tsv',sep='\t',index=False)
    print(pd.DataFrame(rows).to_string(index=False),flush=True);print(delta.to_string(index=False),flush=True)
    dest=ROOT/'results/independent_validation_20260915'
    desktop=__import__('pathlib').Path('<USER_HOME>/Desktop/免疫表达轴/independent_validation_source_data_20260915/gene_disjoint');desktop.mkdir(exist_ok=True)
    for p in out.iterdir():shutil.copy2(p,desktop/p.name)
if __name__=='__main__':main()

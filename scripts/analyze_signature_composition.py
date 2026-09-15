"""Same-patient score structure under compartment and composition controls."""
from pathlib import Path
import json, sys
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tmp/python_packages'))
import gseapy as gp
OUT=ROOT/'results/signature_revision_20260915/composition'
SRC=ROOT/'results/published_signature_external_audit'
CFG=json.loads((ROOT/'config/signature_composition_revision_20260915.json').read_text(encoding='utf-8'))
GROUPS=list(CFG['primary_groups'])
CAT=pd.read_csv(SRC/'signature_catalog.tsv',sep='\t')
SETS={r.signature:r.genes.split(';') for r in CAT.itertuples()}
NAMES=list(SETS)
RNG=np.random.default_rng(CFG['seed'])

def bh(p):
    p=np.asarray(p); order=np.argsort(p); out=np.empty_like(p)
    out[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return out

def structure(x,ref=None):
    r=rankdata(np.asarray(x),axis=0); sd=r.std(axis=0,ddof=1)
    z=(r-r.mean(axis=0))/np.where(sd>0,sd,1)
    _,s,v=np.linalg.svd(z,full_matrices=False); loading=v[0]
    if (loading.sum() if ref is None else loading@ref)<0: loading=-loading
    cor=z.T@z/max(1,len(z)-1); pair=cor[np.triu_indices(x.shape[1],1)]
    return dict(pc1_variance=float(s[0]**2/(s@s)),median_abs_rho=float(np.median(np.abs(pair))),
                effective_dimensions=float((s@s)**2/np.sum(s**4))), z@loading,loading,cor

def score(expr):
    # Constant gene universe and ES, so each score does not depend on other patient profiles.
    frame=pd.DataFrame(np.log1p(expr),index=GENES,columns=CURRENT_IDS)
    run=gp.ssgsea(data=frame.reset_index(names='Gene'),gene_sets=SETS,outdir=None,no_plot=True,
        sample_norm_method='rank',min_size=2,max_size=500,permutation_num=0,threads=1,seed=CFG['seed'],verbose=False)
    return run.res2d.pivot(index='Name',columns='Term',values='ES').astype(float).loc[CURRENT_IDS,NAMES]

def main():
    global GENES,CURRENT_IDS
    counts=pd.read_csv(OUT/'compartment_counts.tsv',sep='\t')
    n=counts.pivot(index='patient_id',columns='compartment',values='n_cells')
    patients=n.index.tolist(); arrays={}
    GENES=np.asarray(sorted(set.intersection(*[set(np.load(OUT/f'counts/{p}.npz')['genes']) for p in patients])))
    gene_audit=[]
    for p in patients:
        d=np.load(OUT/f'counts/{p}.npz'); pos=pd.Index(d['genes']).get_indexer(GENES)
        assert (pos>=0).all()
        arrays[p]={str(g):d['counts'][pos,j] for j,g in enumerate(d['groups'])}
        gene_audit.append(dict(patient_id=p,original_genes=len(d['genes']),common_genes=len(GENES),excluded=';'.join(sorted(set(d['genes'])-set(GENES)))))
    pd.DataFrame(gene_audit).to_csv(OUT/'gene_universe_audit.tsv',sep='\t',index=False)
    original=pd.read_csv(SRC/'multicohort_ssgsea_scores.tsv',sep='\t')
    y=original[original.cohort.eq('GSE288199')].set_index('sample_key').continuous_outcome
    score_rows=[]; coverage=[]; metrics=[]; responses=[]; paired=[]; associations=[]
    all_scores={}
    # No cohort-comparison selection on outcome: all compartment-qualified patients retained.
    for compartment in ['all_QC','T_including_mixed','T_clean','CD8_clean']+GROUPS:
        CURRENT_IDS=n.index[n[compartment]>=50].tolist()
        expr=np.column_stack([arrays[p][compartment]/arrays[p][compartment].sum()*1e6 for p in CURRENT_IDS])
        scored=score(expr); all_scores[compartment]=scored
        for signature,genes in SETS.items():
            present=np.isin(GENES,genes)&(expr.max(axis=1)>0)
            coverage.append(dict(compartment=compartment,signature=signature,n_present=int(present.sum()),n_genes=len(genes),coverage=present.sum()/len(genes)))
        scored.assign(patient_id=scored.index,compartment=compartment,analysis_set='all_eligible_50').to_csv(OUT/f'scores_{compartment}.tsv',sep='\t',index=False)
        score_rows.append(scored.assign(patient_id=scored.index,compartment=compartment,analysis_set='all_eligible_50'))
        m,pc,loading,cor=structure(scored.to_numpy()); metrics.append(dict(analysis_set='all_eligible_50',compartment=compartment,n=len(scored),**m))
        for signature in NAMES:
            test=spearmanr(scored[signature],y.loc[scored.index])
            responses.append(dict(analysis_set='all_eligible_50',compartment=compartment,signature=signature,n=len(scored),rho=test.statistic,p=test.pvalue))
        print(compartment,len(scored),m,flush=True)
    original_scores=original[original.cohort.eq('GSE288199')].set_index('sample_key')[NAMES]
    agreement=[dict(signature=s,rho_original_vs_all_QC=spearmanr(original_scores.loc[patients,s],all_scores['all_QC'].loc[patients,s]).statistic,
        rho_original_vs_clean_T=spearmanr(original_scores.loc[patients,s],all_scores['T_clean'].loc[patients,s]).statistic) for s in NAMES]
    pd.DataFrame(agreement).to_csv(OUT/'score_QC_agreement.tsv',sep='\t',index=False)
    for minimum in CFG['minimum_cells']:
        CURRENT_IDS=n.index[(n[GROUPS]>=minimum).all(axis=1)].tolist()
        ids=CURRENT_IDS; tag=f'matched_four_{minimum}'; count=len(ids)
        if count<8: continue
        raw=np.stack([[arrays[p][g] for g in GROUPS] for p in ids]) # patient, group, gene
        lib=raw.sum(axis=2); weights=lib/lib.sum(axis=1,keepdims=True)
        cpm=raw/lib[:,:,None]*1e6
        native=raw.sum(axis=1)/raw.sum(axis=(1,2))[:,None]*1e6
        fixed=cpm.mean(axis=1)
        ref=(cpm.sum(axis=0)[None,:,:]-cpm)/(count-1)
        composition=(ref*weights[:,:,None]).sum(axis=1)
        assert np.allclose(native,(cpm*weights[:,:,None]).sum(axis=1))
        exprs={'native_four':native.T,'fixed_equal_RNA_four':fixed.T,'composition_only_LOO':composition.T}
        current={c:all_scores[c].loc[ids] for c in ['T_clean','CD8_clean']+GROUPS}
        for name,ex in exprs.items():
            current[name]=score(ex)
            pd.DataFrame(np.log1p(ex),index=GENES,columns=ids).to_csv(OUT/f'{tag}_{name}_logCPM.tsv.gz',sep='\t')
        cellweights=n.loc[ids,GROUPS].div(n.loc[ids,GROUPS].sum(axis=1),axis=0)
        pd.DataFrame(weights,index=ids,columns=GROUPS).rename_axis('patient_id').to_csv(OUT/f'{tag}_RNA_weights.tsv',sep='\t')
        cellweights.rename_axis('patient_id').to_csv(OUT/f'{tag}_cell_weights.tsv',sep='\t')
        native_metrics,_,native_loading,_=structure(current['native_four'].to_numpy())
        fixed_metrics,_,fixed_loading,_=structure(current['fixed_equal_RNA_four'].to_numpy())
        boot=[]
        for _ in range(2000):
            idx=RNG.integers(0,count,count)
            m1,*_=structure(current['native_four'].to_numpy()[idx]); m2,*_=structure(current['fixed_equal_RNA_four'].to_numpy()[idx])
            boot.append([m2[k]-m1[k] for k in ['pc1_variance','median_abs_rho','effective_dimensions']])
        for j,k in enumerate(['pc1_variance','median_abs_rho','effective_dimensions']):
            low,high=np.quantile(np.asarray(boot)[:,j],[.025,.975])
            paired.append(dict(analysis_set=tag,n=count,metric=k,contrast='fixed_equal_RNA_minus_native',difference=fixed_metrics[k]-native_metrics[k],ci_low=low,ci_high=high,n_boot=2000))
        for c,scored in current.items():
            score_rows.append(scored.assign(patient_id=scored.index,compartment=c,analysis_set=tag))
            m,pc,loading,cor=structure(scored.to_numpy()); metrics.append(dict(analysis_set=tag,compartment=c,n=count,**m))
            pd.DataFrame(cor,index=NAMES,columns=NAMES).to_csv(OUT/f'{tag}_{c}_correlation.tsv',sep='\t')
            test=spearmanr(pc,y.loc[ids]); responses.append(dict(analysis_set=tag,compartment=c,signature='PC1',n=count,rho=test.statistic,p=test.pvalue))
            for signature in NAMES:
                test=spearmanr(scored[signature],y.loc[ids])
                responses.append(dict(analysis_set=tag,compartment=c,signature=signature,n=count,rho=test.statistic,p=test.pvalue))
        _,pc,_,_=structure(current['T_clean'].to_numpy())
        for g in GROUPS:
            full_T_RNA=np.array([arrays[p]['T_clean'].sum() for p in ids])
            for kind,w in [('RNA',lib[:,GROUPS.index(g)]/full_T_RNA),('cell',(n.loc[ids,g]/n.loc[ids,'T_clean']).to_numpy())]:
                test=spearmanr(pc,w)
                associations.append(dict(analysis_set=tag,axis='T_clean_PC1',group=g,weight_kind=kind,denominator='all_clean_T',n=count,rho=test.statistic,p=test.pvalue))
        print(tag,count,'native',native_metrics,'fixed',fixed_metrics,flush=True)
    pd.concat(score_rows).to_csv(OUT/'all_compartment_scores.tsv',sep='\t',index=False)
    pd.DataFrame(metrics).to_csv(OUT/'redundancy_metrics.tsv',sep='\t',index=False)
    pd.DataFrame(paired).to_csv(OUT/'paired_structure_changes.tsv',sep='\t',index=False)
    pd.DataFrame(coverage).to_csv(OUT/'signature_coverage.tsv',sep='\t',index=False)
    resp=pd.DataFrame(responses); resp['q_global_exploratory']=bh(resp.p); resp.to_csv(OUT/'response_associations.tsv',sep='\t',index=False)
    assoc=pd.DataFrame(associations); assoc['q_all_composition_tests']=bh(assoc.p); assoc.to_csv(OUT/'axis_composition_associations.tsv',sep='\t',index=False)
    (OUT/'analysis_audit.json').write_text(json.dumps(dict(gseapy=gp.__version__,seed=CFG['seed'],score='ES',gene_universe=len(GENES),
        n_boot=2000,status='exploratory',composition_only='LOO cohort-reference descriptive only; not independent patient predictions',
        bootstrap='paired patient resampling refits rank PCA; fixed equal weights do not require fitting'),indent=2),encoding='utf-8')

if __name__=='__main__': main()

"""Frozen external biological-reference and small clinical-cohort analyses."""
from pathlib import Path
import sys,json,gzip,itertools,hashlib
import numpy as np
import pandas as pd
from scipy.stats import rankdata,spearmanr,pearsonr
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tmp/python_packages'))
import gseapy as gp
from analyze_signature_composition import structure,bh,SETS,NAMES
D=ROOT/'data/external/signature_validation_20260915'
O=ROOT/'results/independent_validation_20260915';O.mkdir(parents=True,exist_ok=True)
CFG=json.loads((ROOT/'config/independent_validation_protocol_20260915.json').read_text())
RNG=np.random.default_rng(CFG['seed'])
def save(d,name): d.to_csv(O/(name+'.tsv'),sep='\t',index=True)
def scores(x,name):
    assert x.index.is_unique and x.columns.is_unique
    coverage=pd.DataFrame([dict(signature=s,present=len(set(g)&set(x.index)),total=len(g),fraction=len(set(g)&set(x.index))/len(g)) for s,g in SETS.items()])
    save(coverage,name+'_coverage');assert coverage.fraction.min()>=.8
    run=gp.ssgsea(data=x.reset_index(names='Gene'),gene_sets=SETS,outdir=None,no_plot=True,sample_norm_method='rank',min_size=2,max_size=500,permutation_num=0,threads=1,seed=CFG['seed'],verbose=False)
    z=run.res2d.pivot(index='Name',columns='Term',values='ES').astype(float).loc[x.columns,NAMES]
    save(z,name+'_scores');return z
def partial(x,y,c):
    a=np.column_stack([np.ones(len(c)),rankdata(c)])
    xr=rankdata(x);yr=rankdata(y)
    return pearsonr(xr-a@np.linalg.lstsq(a,xr,rcond=None)[0],yr-a@np.linalg.lstsq(a,yr,rcond=None)[0]).statistic
def auc(x,y):
    r=rankdata(x);p=int(y.sum());return (r[y==1].sum()-p*(p+1)/2)/(p*(len(y)-p))
def main():
    audit={'protocol_sha256':hashlib.sha256((ROOT/'config/independent_validation_protocol_20260915.json').read_bytes()).hexdigest()}
    clin=pd.read_csv(D/'HS_CPTAC_HNSCC_CLI.tsi',sep='\t',index_col=0).drop(index='CAT')
    rna=pd.read_csv(D/'HS_CPTAC_HNSCC_RNAseq_RSEM_UQ_log2_Tumor.cct',sep='\t',index_col=0)
    prot=pd.read_csv(D/'HS_CPTAC_HNSCC_Proteomics_TMT_Gene_level_Tumor.cct',sep='\t',index_col=0)
    ids=sorted(set(clin.index)&set(rna.columns)&set(prot.columns)); assert len(ids)==108
    audit['CPTAC']={'n':len(ids),'matrix_only_ids':sorted((set(rna.columns)|set(prot.columns))-set(ids)),'RNA_missing':int(rna[ids].isna().sum().sum())}
    # RSEM log2 data already transformed; within-sample rank scoring needs no inversion.
    x=rna[ids].dropna(how='all');assert not x.isna().any().any()
    s=scores(x,'CPTAC');m,pc,loading,cor=structure(s.values)
    save(pd.DataFrame(cor,index=NAMES,columns=NAMES),'CPTAC_signature_correlations')
    union=set().union(*[set(g) for g in SETS.values()])
    refs=[];keep=[]
    for g in CFG['CPTAC']['state_candidates']:
        n=int(prot.loc[g,ids].notna().sum()) if g in prot.index else 0
        use=g not in union and n>=.8*len(ids)
        refs.append(dict(gene=g,n=n,signature_overlap=g in union,included=use))
        if use:keep.append(g)
    assert len(keep)>=3
    save(pd.DataFrame(refs),'CPTAC_reference_gene_audit')
    p=prot.loc[keep,ids].T
    state=((p-p.mean())/p.std(ddof=1)).mean(axis=1,skipna=False)
    table=pd.DataFrame({'PC1':pc,'CD3_IHC':pd.to_numeric(clin.loc[ids,'CD3_IHC_count'],errors='coerce'),'PTPRC_protein':prot.loc['PTPRC',ids],'state_protein':state},index=ids)
    save(table,'CPTAC_patient_references');rows=[]
    for label,ref,adjust in [('CD3_IHC','CD3_IHC',False),('PTPRC_protein','PTPRC_protein',False),('state_protein_given_CD3','state_protein',True)]:
        valid=table[['PC1',ref]+(['CD3_IHC'] if adjust else [])].notna().all(axis=1)
        j=np.flatnonzero(valid);t=table.iloc[j];v=partial(t.PC1,t[ref],t.CD3_IHC) if adjust else spearmanr(t.PC1,t[ref]).statistic
        if adjust:
            from scipy.stats import t as tdist
            pv=2*tdist.sf(abs(v)*np.sqrt((len(j)-3)/(1-v*v)),len(j)-3)
        else:pv=spearmanr(t.PC1,t[ref]).pvalue
        boots=[]
        # Resample all 108; refit axis and then select observed reference complete cases.
        for _ in range(CFG['bootstrap']):
            ix=RNG.integers(0,len(ids),len(ids));_,bpc,_,_=structure(s.values[ix],ref=loading)
            bt=table.iloc[ix].copy();bt['PC1']=bpc
            if adjust:
                bp=p.iloc[ix];bt['state_protein']=((bp-bp.mean())/bp.std(ddof=1)).mean(axis=1,skipna=False).to_numpy()
            bt=bt.dropna(subset=['PC1',ref]+(['CD3_IHC'] if adjust else []))
            boots.append(partial(bt.PC1,bt[ref],bt.CD3_IHC) if adjust else spearmanr(bt.PC1,bt[ref]).statistic)
        lo,hi=np.quantile(boots,[.025,.975]);rows.append(dict(test=label,n=len(t),effect=v,ci_low=lo,ci_high=hi,p=pv))
    refs=pd.DataFrame(rows);refs['q']=bh(refs.p);save(refs,'CPTAC_associations');audit['CPTAC'].update(m)
    print(refs.to_string(index=False),flush=True)
    soft=gzip.open(D/'GSE284162_family.soft.gz','rt',encoding='utf-8').read();records=[]
    for block in soft.split('^SAMPLE = ')[1:]:
        rec={'gsm':block.splitlines()[0]}
        for line in block.splitlines():
            if line.startswith('!Sample_characteristics_ch1 = '):
                k,v=line.split(' = ',1)[1].split(': ',1);rec[k]=v
            if line.startswith('!Sample_platform_id = '):rec['platform']=line.split(' = ')[1]
        records.append(rec)
    meta=pd.DataFrame(records).set_index('individual').sort_index();save(meta,'GSE284162_metadata')
    assert meta.index.is_unique and meta.ID.is_unique
    counts=pd.read_excel(D/'GSE284162_raw_counts.xlsx',index_col=0)
    audit['GSE284162']={'original_genes':len(counts),'duplicate_gene_rows':int(counts.index.duplicated().sum())}
    assert counts.index.notna().all() and not counts.isna().any().any() and (counts.values>=0).all()
    counts=counts.groupby(level=0).sum().loc[:,meta.index];counts=counts.loc[counts.sum(axis=1)>0]
    y=meta.mpr.eq('Yes').astype(int).to_numpy();assert y.sum()==4 and len(y)==12
    save(pd.crosstab(meta.platform,meta.mpr),'GSE284162_platform_response')
    s=scores(np.log1p(counts/counts.sum()*1e6),'GSE284162');m,pc,loading,cor=structure(s.values)
    audit['GSE284162'].update(m);save(pd.DataFrame(cor,index=NAMES,columns=NAMES),'GSE284162_signature_correlations')
    s=s.assign(PC1=pc);save(s.join(meta),'GSE284162_patient_scores_metadata')
    permutations=[]
    for comb in itertools.combinations(range(12),4):
        z=np.zeros(12,dtype=int);z[list(comb)]=1;permutations.append(z)
    restricted=[z for z in permutations if all(z[meta.platform.eq(b)].sum()==y[meta.platform.eq(b)].sum() for b in meta.platform.unique())]
    audit['GSE284162']['permutations']={'all':len(permutations),'platform_restricted':len(restricted)}
    rows=[]
    for name in CFG['GSE284162']['tests']:
        v=auc(s[name].values,y);stats={key:np.mean([abs(auc(s[name].values,z)-.5)>=abs(v-.5)-1e-12 for z in perms]) for key,perms in [('p_exact',permutations),('p_platform',restricted)]}
        boots=[]
        for _ in range(CFG['bootstrap']):
            ix=np.concatenate([RNG.choice(np.flatnonzero(y==a),sum(y==a),replace=True) for a in [0,1]])
            bx=structure(s[NAMES].values[ix],ref=loading)[1] if name=='PC1' else s[name].values[ix]
            boots.append(auc(bx,y[ix]))
        lo,hi=np.quantile(boots,[.025,.975]);rows.append(dict(test=name,n=12,n_positive=4,AUC=v,ci_low=lo,ci_high=hi,**stats))
    tests=pd.DataFrame(rows)
    for field in ['p_exact','p_platform']:tests[field.replace('p_','q_')]=bh(tests[field])
    save(tests,'GSE284162_response_associations');print(tests.to_string(index=False),flush=True)
    (O/'analysis_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8');print(json.dumps(audit,indent=2),flush=True)
if __name__=='__main__':main()

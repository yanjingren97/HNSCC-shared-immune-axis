from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import rankdata
R=Path(__file__).resolve().parents[1]
O=R/'results/biological_specificity_20260917';O.mkdir(exist_ok=True)
D=Path.home()/'Desktop/免疫表达轴/生物学特异性修订_20260917';D.mkdir(exist_ok=True)
CFG=json.loads((R/'config/biological_specificity_20260917.json').read_text())
CAT=pd.read_csv(R/'results/published_signature_external_audit/signature_catalog.tsv',sep='\t')
SETS={r.signature:r.genes.split(';') for r in CAT.itertuples()};NAMES=list(SETS)
ALIASES=dict(zip(NAMES,['GZMK40','CYT2','EXH5','IFNG6','TLS9','GEP18','Korea30','Dysf238','Shared84']))
SHARED=NAMES[-1]
def standardized(x):
    x=rankdata(np.asarray(x),axis=0);sd=x.std(axis=0,ddof=1)
    return (x-x.mean(axis=0))/np.where(sd>0,sd,1)
def pc(x):
    z=standardized(x);_,s,v=np.linalg.svd(z,full_matrices=False);v=v[0]
    if v.sum()<0:v=-v
    return float(s[0]**2/(s@s)),z@v
def load_expr(name):
    if name=='CPTAC':
        ids=pd.read_csv(R/'results/independent_validation_20260915/CPTAC_patient_references.tsv',sep='\t',index_col=0).index
        return pd.read_csv(R/'data/external/signature_validation_20260915/HS_CPTAC_HNSCC_RNAseq_RSEM_UQ_log2_Tumor.cct',sep='\t',index_col=0)[ids]
    comp='native_four' if 'native' in name else 'fixed_equal_RNA_four'
    return pd.read_csv(R/f'results/signature_revision_20260915/composition/matched_four_50_{comp}_logCPM.tsv.gz',sep='\t',index_col=0)
class FastSSGSEA:
    """Closed form of integrated ssGSEA running sum with rank^.25 weights."""
    def __init__(self,x):
        self.genes=x.index;self.lookup={g:i for i,g in enumerate(x.index)};self.n=len(x)
        ranks=rankdata(x.values,axis=0);self.weight=ranks**.25
        # Reuse the installed GSEApy engine's exact ordering of tied genes.
        # A different tie order notably changes scores for sparse transcripts.
        import sys
        sys.path.insert(0,str(R/'tmp/python_packages'))
        from gseapy.gse import ssgsea_rs,CorrelType
        dummy=next(g for g in SETS.values() if len(self.ids(g))>=2)
        raw=ssgsea_rs(x.index.tolist(),(ranks*10000/self.n).tolist(),{'order_only':dummy},.25,2,500,0,CorrelType.Rank,1,CFG['seed'])
        order=np.asarray(raw.indices).T
        pos=np.empty_like(order);np.put_along_axis(pos,order,np.arange(self.n)[:,None],axis=0)
        self.integral=self.n-pos
    def ids(self,gs):return np.array([self.lookup[g] for g in gs if g in self.lookup],dtype=int)
    def score_ids(self,idx):
        k=len(idx);w=self.weight[idx];p=self.integral[idx]
        return (w*p).sum(axis=0)/w.sum(axis=0)-(self.n*(self.n+1)/2-p.sum(axis=0))/(self.n-k)
    def score(self,sets):return np.column_stack([self.score_ids(self.ids(g)) for g in sets.values()])

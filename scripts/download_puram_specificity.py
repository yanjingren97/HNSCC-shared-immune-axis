import urllib.request,gzip
from pathlib import Path
p=Path('data/external/specificity_20260917/GSE103322_HNSCC_all_data.txt.gz')
if not p.exists():urllib.request.urlretrieve('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103322/suppl/'+p.name,p)
with gzip.open(p,'rt') as f:
 for i in range(9):
  a=f.readline().rstrip().split('\t');print(i,len(a),a[:12])

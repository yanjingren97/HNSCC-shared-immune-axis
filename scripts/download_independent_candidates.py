from pathlib import Path
import re,urllib.request,concurrent.futures
D=Path('data/external/ifng6_stability_20260918');s=(D/'engraftment.html').read_text(encoding='utf-8')
for m in re.finditer('MOESM[34567]_ESM.xlsx',s):print(re.sub('<[^>]+>',' ',s[max(0,m.start()-100):m.end()+200])[-350:])
urls=sorted(set(re.findall(r'https[^"<> ]+MOESM[34567]_ESM.xlsx',s)))
urls+=['https://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/GSE288406/suppl/GSE288406_HNSCC_IPA_normalised_lcpm.csv.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/GSE288406/suppl/GSE288406_HNSCC_IPA_normalised_metadata.csv.gz']
def fetch(u):
 p=D/u.split('/')[-1]
 if not p.exists():urllib.request.urlretrieve(u,p)
 return p.name,p.stat().st_size
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:print(list(pool.map(fetch,urls)))

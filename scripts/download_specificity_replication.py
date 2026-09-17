from pathlib import Path
import urllib.request,concurrent.futures,json,hashlib,re
D=Path('data/external/specificity_20260917')
def get(url):
 with urllib.request.urlopen(url,timeout=180) as r:return r.read()
p=get('https://gdc.cancer.gov/about-data/publications/panimmune').decode()
url=next(u for u,t in re.findall(r'<a[^>]+href="([^"]+)"[^>]*>(.*?)</a>',p,re.S) if 'TCGA_all_leuk' in t)
items={'TCGA_leukocyte_methylation.tsv':url,'TCGA_HNSC_HiSeqV2.gz':'https://tcga.xenahubs.net/download/TCGA.HNSC.sampleMap/HiSeqV2.gz','TCGA_HNSC_HiSeqV2.json':'https://tcga.xenahubs.net/download/TCGA.HNSC.sampleMap/HiSeqV2.json'}
def fetch(kv):
 name,url=kv;b=get(url);(D/name).write_bytes(b);return dict(file=name,url=url,bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
 records=list(pool.map(fetch,items.items()))
(D/'replication_downloads.json').write_text(json.dumps(records,indent=2));print(records)

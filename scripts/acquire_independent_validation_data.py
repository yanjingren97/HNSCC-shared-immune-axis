from pathlib import Path
import urllib.request,urllib.parse,hashlib,json,gzip
from html.parser import HTMLParser
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/external/signature_validation_20260915'
class Links(HTMLParser):
    def __init__(self): super().__init__(); self.links=[]
    def handle_starttag(self,tag,attrs):
        if tag=='a': self.links.extend(v for k,v in attrs if k=='href')
def main():
    OUT.mkdir(parents=True,exist_ok=True)
    page='https://www.linkedomics.org/data_download/CPTAC-HNSCC/'
    def fetch(url):
        return urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'}),timeout=60).read()
    try: html=fetch(page)
    except Exception:
        page='https://linkedomics.org/data_download/CPTAC-HNSCC/'
        html=fetch(page)
    (OUT/'CPTAC_download_page.html').write_bytes(html)
    parser=Links();parser.feed(html.decode())
    links=[urllib.parse.urljoin(page,x) for x in parser.links]
    print('\n'.join(links),flush=True)
    selected=[x for x in links if any(k in x.lower() for k in ['clinical','rna','proteome','_cli','_prot'])]
    (OUT/'available_links.json').write_text(json.dumps(selected,indent=2),encoding='utf-8')
    downloads={
      'GSE284162_family.soft.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE284nnn/GSE284162/soft/GSE284162_family.soft.gz',
      'GSE284162_raw_counts.xlsx':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE284nnn/GSE284162/suppl/GSE284162_raw_counts.xlsx'}
    for x in selected:
        name=urllib.parse.unquote(x.rsplit('/',1)[-1])
        low=name.lower()
        if 'clinical' in low or '_cli' in low or ('tumor' in low and ('proteome' in low or ('rna' in low and 'isoform' not in low and 'mirna' not in low and 'circrna' not in low))):
            downloads[name]=x
    records=[]
    for name,url in downloads.items():
        p=OUT/name
        try:
            if not p.exists():
                p.write_bytes(fetch(url))
            records.append(dict(file=name,url=url,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),status='downloaded'))
            print(name,p.stat().st_size,flush=True)
        except Exception as ex:
            records.append(dict(file=name,url=url,status='failed',error=str(ex)));print(name,str(ex),flush=True)
    (OUT/'download_manifest.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
if __name__=='__main__':main()

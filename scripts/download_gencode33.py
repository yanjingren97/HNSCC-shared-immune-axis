from pathlib import Path
import urllib.request,gzip,re,json
D=Path('data/external/ifng6_stability_20260918');p=D/'gencode.v33.annotation.gtf.gz'
if not p.exists():
 with urllib.request.urlopen('https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gtf.gz',timeout=45) as r,p.open('wb') as f:
  while chunk:=r.read(1024*1024):f.write(chunk)
rows=[]
with gzip.open(p,'rt') as f:
 for line in f:
  if line.startswith('#'):continue
  col=line.rstrip().split('\t')
  if col[2]=='gene':
   a=dict(re.findall(r'(\w+) "([^"]+)"',col[8]));rows.append((a['gene_id'].split('.')[0],a['gene_name']))
(D/'gencode_v33_gene_symbols.tsv').write_text('ENSG\tgene\n'+'\n'.join('\t'.join(r) for r in rows));print(len(rows))

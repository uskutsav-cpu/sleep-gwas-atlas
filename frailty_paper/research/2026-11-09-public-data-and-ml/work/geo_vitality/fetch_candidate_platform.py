import requests,gzip,io,csv,json
from pathlib import Path
u='https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL15nnn/GPL15988/soft/GPL15988_family.soft.gz';r=requests.get(u,stream=True,timeout=45);r.raise_for_status();targets={'FOXP2','SMAD5','RBM6','MEIS1'};rows=[]
with gzip.GzipFile(fileobj=r.raw) as z:
 f=io.TextIOWrapper(z)
 for line in f:
  if line.startswith('!platform_table_begin'):break
 header=next(f).rstrip('\n').split('\t');assert header[:4]==['ID','GB_ACC','GENE SYMBOL','GENE DESCRIPTION'];n=0
 for line in f:
  if line.startswith('!platform_table_end'):break
  x=line.rstrip('\n').split('\t');n+=1
  if len(x)>=4 and x[2] in targets:rows.append({'probe':x[0],'accession':x[1],'gene':x[2],'description':x[3]})
Path('outputs/stage6_gse40366_candidate_probe_map.json').write_text(json.dumps({'platform':'GPL15988','source':u,'rows':rows,'n_platform_rows':n,'scope':'Exact submitted gene symbols; annotation may be outdated or cross-hybridize'},indent=2)+'\n');print('platform rows',n,'candidate probes',rows)

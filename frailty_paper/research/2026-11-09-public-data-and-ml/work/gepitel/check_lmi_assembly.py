import concurrent.futures,json,requests
from pathlib import Path
cases=[('rs34290943','3',49161660),('rs8034223','15',27044693),('rs182813191','9',14419143),('rs9638291','7',70690615)]
def fetch(case):
 rs,chrom,pos=case;out={'variant':rs,'workbook_chr':chrom,'workbook_pos':pos,'mappings':{}}
 for key,host in [('GRCh37','grch37.rest.ensembl.org'),('GRCh38','rest.ensembl.org')]:
  url='https://'+host+'/variation/human/'+rs;r=requests.get(url,params={'content-type':'application/json'},timeout=30);r.raise_for_status();out['mappings'][key]=[{'chr':x['seq_region_name'],'position':x['start'],'alleles':x['allele_string']} for x in r.json()['mappings'] if x.get('coord_system')=='chromosome' and x.get('seq_region_name')==chrom]
 return out
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:data=list(pool.map(fetch,cases))
Path('outputs/stage6_gepitel_lmi_build_check.json').write_text(json.dumps(data,indent=2)+'\n');print('checked',len(data),'37 matches',sum(any(m['position']==x['workbook_pos'] for m in x['mappings']['GRCh37']) for x in data),'38 matches',sum(any(m['position']==x['workbook_pos'] for m in x['mappings']['GRCh38']) for x in data))

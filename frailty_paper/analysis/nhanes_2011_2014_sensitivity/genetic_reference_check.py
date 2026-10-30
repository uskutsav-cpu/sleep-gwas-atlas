#!/usr/bin/env python3
"""GRCh37 reference inventory, not target-genotype or GWAS effect verification."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,time
import pandas as pd,requests
P=Path(__file__).resolve().parent;cache=P/'reference_grch37';cache.mkdir(exist_ok=True)
scores=pd.concat([pd.read_csv(P/f'score_{n}_variants_standardized.csv').assign(score=n) for n in [15,22]],ignore_index=True)
def get(rs):
 path=cache/f'{rs}.json';url=f'https://grch37.rest.ensembl.org/variation/human/{rs}?content-type=application/json'
 if path.exists():return rs,url,200,path.read_bytes()
 for attempt in range(3):
  try:
   r=requests.get(url,timeout=25)
   if r.status_code==200:
    path.write_bytes(r.content);return rs,url,r.status_code,r.content
  except requests.RequestException as e:status=str(e)
  else:status=f'HTTP {r.status_code}'
  time.sleep(1+attempt)
 return rs,url,status,b''
with ThreadPoolExecutor(max_workers=3) as pool:response=dict((r[0],r[1:]) for r in pool.map(get, sorted(scores.rsid.unique())))
rows=[]
for _,r in scores.iterrows():
 url,status,body=response[r.rsid];result={'score':r.score,'rsid':r.rsid,'chromosome_expected':str(r.chromosome),'position_expected':int(r.position),'counted_allele':r.effect_counted_allele,'other_allele':r.other_allele,'source_GWAS_build':r.genome_build_reported,'source_build_status':r.build_certainty,'reference_url':url,'response_sha256':hashlib.sha256(body).hexdigest() if body else '', 'request_status':status}
 if body:
  data=json.loads(body);maps=[m for m in data.get('mappings',[]) if m.get('assembly_name')=='GRCh37' and str(m.get('seq_region_name'))==str(r.chromosome) and m.get('start')==int(r.position) and m.get('end')==int(r.position)];result['exact_locus_mappings']=len(maps);als=set(maps[0].get('allele_string','').split('/')) if maps else set();result['reference_alleles']='/'.join(sorted(als));result['score_pair_present']=bool({r.effect_counted_allele,r.other_allele}<=als);result['reference_multiallelic']=len(als)>2;result['reference_build']='GRCh37'
 else:result.update(exact_locus_mappings=None,reference_alleles='',score_pair_present=None,reference_multiallelic=None,reference_build='unknown')
 result['NHANES_2011_2012_genotype_coverage']='unknown: CDC catalog does not list this cycle';result['FFS_weight_direction']='cached A1 alignment inferred; target and source contrast unresolved';rows.append(result)
out=pd.DataFrame(rows);out.to_csv(P/'supplement_genetic_reference_audit.csv',index=False)
print('score rows',len(out),'requests ok',sum(out.request_status==200),'coordinate matches',sum(out.exact_locus_mappings==1),'pairs present',sum(out.score_pair_present==True),'failed',sum(out.request_status!=200));print(out[out.request_status!=200][['score','rsid','request_status']].to_string(index=False))

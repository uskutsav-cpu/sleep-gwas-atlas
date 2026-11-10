#!/usr/bin/env python3
import csv, hashlib,json,subprocess,concurrent.futures
from pathlib import Path
previous=list(csv.DictReader(open('outputs/ncds_22_grch38_ensembl_audit.csv')))
failed={r['variant'] for r in previous if r['fetch_status']!='ok'}
rows=[r for r in csv.DictReader(open('outputs/ncds_score_variant_request.csv')) if r['score']=='22_count' and r['variant'] in failed]
def fetch(r):
 url='https://rest.ensembl.org/variation/human/'+r['variant']+'?content-type=application/json'
 cp=subprocess.run(['curl','-fsSL','--max-time','10',url],capture_output=True)
 result=dict(variant=r['variant'],chr38=r['chr38'],pos38=r['pos38'],fetch_status='error',matching_mappings='',reference_alleles='',both_score_alleles_listed='',multiallelic='',url=url,response_sha256='')
 if cp.returncode: result['fetch_status']='curl_exit_'+str(cp.returncode); return result
 try: x=json.loads(cp.stdout)
 except Exception:result['fetch_status']='invalid_json'; return result
 hit=[v for v in x.get('mappings',[]) if v.get('assembly_name')=='GRCh38' and str(v.get('seq_region_name'))==r['chr38'] and str(v.get('start'))==r['pos38'] and str(v.get('end'))==r['pos38']]
 result.update(fetch_status='ok',matching_mappings=len(hit),response_sha256=hashlib.sha256(cp.stdout).hexdigest())
 if len(hit)==1:
  alleles=hit[0].get('allele_string',''); aset=set(alleles.split('/'))
  result.update(reference_alleles=alleles,both_score_alleles_listed={r['counted_allele_38_candidate'],r['other_allele_38_candidate']}<=aset,multiallelic=len(aset)>2)
 return result
with concurrent.futures.ThreadPoolExecutor(max_workers=12) as p: results=list(p.map(fetch,rows))
out=Path('outputs/ncds_22_grch38_ensembl_retry.csv')
with out.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=results[0]);w.writeheader();w.writerows(results)
print(json.dumps({'requested':len(results),'retrieved':sum(x['fetch_status']=='ok' for x in results),'matching_locus':sum(x['matching_mappings']==1 for x in results),'pair_listed':sum(x['both_score_alleles_listed'] is True for x in results),'multiallelic':sum(x['multiallelic'] is True for x in results),'failed':[x['variant'] for x in results if x['fetch_status']!='ok']},indent=2))

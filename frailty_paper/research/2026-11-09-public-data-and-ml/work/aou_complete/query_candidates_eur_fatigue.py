import csv,json,subprocess,time,hashlib,urllib.parse
from pathlib import Path
p=Path('outputs/stage6_kunorozva_meta_candidate_extract.csv')
with p.open(newline='') as f:candidates=list(csv.DictReader(f))
assert len(candidates)==33 and len({r['candidate'] for r in candidates})==33
out=Path('outputs/stage6_aou_eur_fatigue_candidate_association.csv')
fields=['candidate','phenotype','chrom','position','ref','alt','status','pvalue','beta','se','af','ancestry','sequencing_type','storage_source','returned_count','error']
prior={}
if out.exists():
 with out.open(newline='') as f:
  for r in csv.DictReader(f):
   if r['status'] in ('matched','no_matching_row'):prior[(r['candidate'],r['phenotype'])]=r
for c in candidates:
 chrom=c['#chrom'];pos=int(c['pos']);ref=c['ref'];alt=c['alt'];vid=f'chr{chrom}-{pos}-{ref}-{alt}'
 for phenotype in ('SS_807',):
  key=(c['candidate'],phenotype)
  if key in prior:continue
  query=urllib.parse.urlencode({'analysis_id':phenotype,'ancestry_group':'eur','sequencing_type':'genome','query_mode':'slow'})
  url=f'https://allbyall.researchallofus.org/api/variants/associations/interval/chr{chrom}:{pos-1}-{pos+1}?{query}'
  proc=subprocess.run(['curl','-fsSL','--max-time','25',url],capture_output=True,text=True)
  r=dict(candidate=c['candidate'],phenotype=phenotype,chrom=chrom,position=pos,ref=ref,alt=alt,status='error',pvalue='',beta='',se='',af='',ancestry='',sequencing_type='',storage_source='',returned_count='',error='')
  if proc.returncode:
   r['error']=proc.stderr.strip()[:200]
  else:
   try:
    j=json.loads(proc.stdout);assert j.get('storage_source')=='hail_gcs' and isinstance(j.get('data'),list)
    matches=[x for x in j['data'] if x['variant_id']==vid]
    assert len(matches)<=1
    r['returned_count']=j['count'];r['storage_source']=j['storage_source']
    if not matches:r['status']='no_matching_row'
    else:
     x=matches[0];assert x['phenotype']==phenotype and x['ancestry']=='eur' and x['sequencing_type']=='genome'
     assert 0<=x['pvalue']<=1 and x['se']>0
     r.update(status='matched',pvalue=x['pvalue'],beta=x['beta'],se=x['se'],af=x['af'],ancestry=x['ancestry'],sequencing_type=x['sequencing_type'])
   except Exception as e:r['error']=repr(e)[:200]
  prior[key]=r
  with out.open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(prior.values())
  print(c['candidate'],phenotype,r['status'],flush=True)
  time.sleep(.15)
summary={'source_file_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'requested_pairs':33,'recorded':len(prior),'status_counts':{s:sum(x['status']==s for x in prior.values()) for s in ('matched','no_matching_row','error')},'metadata':'CDR v8; All by All browser 0.1.0-317fe0c-202609021906; genome, EUR genetic-similarity group, slow Hail interval','caveat':'Only summary-level; association data do not measure performance of project score in individuals. Empty or error is not a null association; inspect exact allele semantics and overlap.'}
Path('outputs/stage6_aou_eur_fatigue_candidate_association_receipt.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)

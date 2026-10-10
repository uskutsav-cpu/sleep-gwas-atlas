#!/usr/bin/env python3
"""Bounded official policy/record metadata GETs only, never source data bodies."""
import datetime,hashlib,json,re,subprocess,sys
from pathlib import Path
P=Path(__file__).resolve().parents[1]
PLAN=P/'source_provenance/source_policy_followup_microfetch_plan_v4_2.json'
RECEIPT=P/'source_provenance/source_policy_followup_microfetch_receipt_v4_2.json'
EXPECTED='f27fa1daa1a7fd198966a44ab5a48bf54b7ca50c41701e4f05d5bfdd3996793e'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def regular(p):
 p=Path(p);assert not p.is_symlink() and all(not q.is_symlink() for q in p.parents);return p
assert sha(regular(PLAN))==EXPECTED
plan=json.loads(PLAN.read_text());cache=regular(plan['cache_namespace'])
assert not cache.exists() and not RECEIPT.exists()
version=subprocess.run(['/usr/bin/curl','-q','--version'],capture_output=True,text=True,check=True).stdout
v=re.search(r'^curl (\d+)\.(\d+)\.(\d+)',version);assert v and tuple(map(int,v.groups()))>=(8,4,0)
cache.mkdir(parents=True);total=0;responses=[]
for i,item in enumerate(plan['requests'],1):
 assert re.fullmatch(r'[A-Za-z0-9_]+',item['id'])
 assert item['url'].startswith('https://') and not re.search(r'\.(gz|bgz|zip|vcf|bed|bim|fam)(\?|$)',item['url'])
 p=regular(cache/(item['id']+'.response.txt'))
 with p.open('xb'):pass
 cmd=['/usr/bin/curl','-q','--location','--max-redirs','3','--proto','=https','--proto-redir','=https','--tlsv1.2','--max-time','30','--connect-timeout','10','--max-filesize',str(plan['per_response_body_cap_bytes']),'--output',str(p),'--silent','--show-error','--write-out','%{http_code}\t%{size_download}\t%{url_effective}\t%{content_type}\n',item['url']]
 run=subprocess.run(cmd,capture_output=True,text=True,timeout=35)
 assert len(run.stdout.encode())<=plan['per_response_header_cap_bytes'] and len(run.stderr.encode())<=plan['per_response_header_cap_bytes']
 size=p.stat().st_size;assert size<=plan['per_response_body_cap_bytes'];total+=size
 assert total<=plan['total_retained_cap_bytes']
 fields=run.stdout.strip().split('\t',3)
 response={**item,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'returncode':run.returncode,'http_code':fields[0] if fields else None,'curl_writeout':run.stdout.strip(),'stderr':run.stderr[:4000],'path':str(p),'bytes':size,'sha256':sha(p),'metadata_only':True,'GWAS_or_reference_body':False}
 responses.append(response)
 if i%5==0 or i==len(plan['requests']):print(json.dumps({'metadata_requests_completed':i,'total':len(plan['requests']),'retained_bytes':total}),flush=True)
receipt={'schema':'bounded_public_source_policy_followup_microfetch_receipt_v4_2','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'plan_sha256':EXPECTED,'script_sha256':sha(__file__),'curl_sha256':sha('/usr/bin/curl'),'curl_version':version,'actual_response_body_bytes_retained':total,'maximum_response_bytes':max(r['bytes'] for r in responses),'request_count':len(responses),'responses':responses,'TLS_verification_not_disabled':True,'automatic_retry':False,'GWAS_or_reference_body_bytes':0,'form_submission_or_investigator_contact':False,'copyrighted_fulltext_committed':False}
with RECEIPT.open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
print(json.dumps({'receipt':str(RECEIPT),'receipt_sha256':sha(RECEIPT),'retained_bytes':total}),flush=True)

#!/usr/bin/env python3
"""Public, bounded acquisition with pre-acquisition manifests and complete pagination.

No summary-statistic GWAS archive, access form, authenticated service or browser bypass.
Raw public literature is retained in ignored work/. Derived metadata is versioned.
"""
import argparse,csv,hashlib,json,re,subprocess,time,urllib.parse
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
AREA=ROOT/'brain6/translational_psychiatry_research_v1/research_v1'
STORE=ROOT/'work/tp-literature-20261007'
LIMIT=30_000_000
CUTOFF='2026-10-07'
QUERY={
 'sleep_psychiatric_genetics':'TITLE_ABS:(sleep OR insomnia OR chronotype) AND TITLE_ABS:(psychiatric OR depression OR schizophrenia OR bipolar OR ADHD) AND TITLE_ABS:(GWAS OR "genome-wide association" OR pleiotropic OR colocalization OR "genetic correlation")',
 'sleep_local_finemapping':'TITLE_ABS:(sleep OR insomnia) AND TITLE_ABS:(LAVA OR SUPERGNOVA OR "local genetic" OR "fine mapping" OR finemapping OR "credible set")',
 'sleep_brain_molecular':'TITLE_ABS:(sleep OR insomnia) AND TITLE_ABS:(eQTL OR sQTL OR "cell type" OR "cell-type" OR "single cell" OR "gene regulation") AND TITLE_ABS:(GWAS OR genetic OR colocalization)',
 'psychiatric_pleiotropy':'TITLE_ABS:(ADHD OR depression OR insomnia) AND TITLE_ABS:(pleiotropic OR colocalization OR "fine mapping" OR "shared effector")',
 'sleep_independent_sources':'TITLE_ABS:(insomnia OR sleep) AND TITLE_ABS:(GWAS OR "genome-wide association") AND TITLE_ABS:(FinnGen OR HUNT OR STARRS OR Partners OR Vanderbilt OR "Million Veteran")',
}

def tsv(path,rows,fields=None):
 path.parent.mkdir(parents=True,exist_ok=True)
 if fields is None:fields=list(rows[0]) if rows else []
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)

class Acquisition:
 def __init__(self,mode,retry=False):
  STORE.mkdir(parents=True,exist_ok=True)
  self.manifest=[];self.receipts=[];self.mode=mode;self.retry=retry
  receipt_path=AREA/f'novelty/{mode}_acquisition_receipts.json'
  self.previous={r['id']:r for r in json.loads(receipt_path.read_text())} if receipt_path.exists() else {}
 def fetch(self,id,url,suffix='.json',purpose='Public metadata',max_bytes=LIMIT):
  path=STORE/(id+suffix)
  row={'id':id,'url':url,'max_bytes':max_bytes,'purpose':purpose,'access':'UNAUTHENTICATED_PUBLIC_SOURCE','planned_path':str(path.relative_to(ROOT))}
  self.manifest.append(row);tsv(AREA/f'novelty/{self.mode}_acquisition_manifest.tsv',self.manifest)
  receipt=dict(row,retrieved_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
  if path.exists():
   payload=path.read_bytes()
   receipt=dict(self.previous.get(id,receipt),**row)
   receipt.update(status='REUSED_LOCAL_ACQUIRED_OBJECT',bytes=len(payload),sha256=hashlib.sha256(payload).hexdigest(),rechecked_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
  elif self.previous.get(id,{}).get('status')=='UNAVAILABLE' and not self.retry:
   receipt=dict(self.previous[id]);payload=None
  else:
   temp=path.with_suffix('.download')
   p=subprocess.run(['curl','--location','--silent','--show-error','--max-time','40','--max-filesize',str(max_bytes),'--output',str(temp),'--write-out','%{json}',url],capture_output=True,text=True)
   try:meta=json.loads(p.stdout)
   except ValueError:meta={}
   receipt.update(http_status=meta.get('http_code'),content_type=meta.get('content_type'),final_url=meta.get('url_effective'),curl_exit=p.returncode)
   payload=temp.read_bytes() if temp.exists() else b''
   if p.returncode or meta.get('http_code')!=200 or len(payload)>max_bytes:
    receipt.update(status='UNAVAILABLE',bytes=len(payload),error=p.stderr.strip());payload=None
    if temp.exists():temp.unlink()
   else:
    temp.rename(path);receipt.update(status='RETRIEVED',bytes=len(payload),sha256=hashlib.sha256(payload).hexdigest())
  self.receipts.append(receipt)
  (AREA/f'novelty/{self.mode}_acquisition_receipts.json').write_text(json.dumps(self.receipts,indent=2)+'\n')
  return payload

def search(a):
 papers={};query_summary=[]
 for name,q in QUERY.items():
  start='2018-01-01' if name=='sleep_independent_sources' else '2025-01-01'
  q+=f' AND FIRST_PDATE:[{start} TO {CUTOFF}]'
  cursor='*';page=0;ids=set();hitcount=None;complete=False
  while True:
   page+=1
   url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode({'query':q,'format':'json','resultType':'core','pageSize':250,'cursorMark':cursor})
   payload=a.fetch(f'{name}_page{page:03}',url,purpose='Focused date-bounded literature metadata; cursor pagination until source reports completion')
   if not payload:break
   d=json.loads(payload);hitcount=d['hitCount'];rr=d.get('resultList',{}).get('result',[])
   for r in rr:
    key=(r.get('source',''),r.get('id',''));ids.add(key)
    papers.setdefault(key,dict(r,queries=[]))['queries'].append(name)
   nextcursor=d.get('nextCursorMark')
   if len(ids)>=hitcount or not rr or nextcursor==cursor:complete=len(ids)==hitcount;break
   if not nextcursor:break
   cursor=nextcursor
  query_summary.append({'query_id':name,'query':q,'pages_attempted':page,'source_hit_count':hitcount,'unique_records_retrieved':len(ids),'complete_for_query':complete,'cutoff':CUTOFF,'interpretation':'Complete indexed query retrieval does not prove exhaustive worldwide literature or novelty'})
  print(name,hitcount,len(ids),complete,flush=True)
 rows=[]
 for key,r in sorted(papers.items()):
  txt=(r.get('title','')+' '+r.get('abstractText','')).lower()
  close=bool(re.search(r'insomnia',txt) and re.search(r'adhd|attention.deficit|depression|depressive',txt) and re.search(r'pleiotrop|colocali|local genetic|fine.?map|shared effector',txt))
  rows.append({'source':key[0],'id':key[1],'doi':r.get('doi',''),'pmid':r.get('pmid',''),'pmcid':r.get('pmcid',''),'title':r.get('title',''),'first_publication_date':r.get('firstPublicationDate',''),'query_ids':';'.join(sorted(set(r['queries']))),'screen':'PRIORITIZE_FULLTEXT_COMPETITOR' if close else 'METADATA_RETRIEVED_NO_AUTOMATIC_NOVELTY_INFERENCE','metadata_review':'TITLE_ABSTRACT_RULE_TRIAGE_NOT_FULLTEXT_REVIEW'})
 tsv(AREA/'novelty/search_records.tsv',rows)
 tsv(AREA/'novelty/query_completion.tsv',query_summary)
 (STORE/'focused_complete_metadata.json').write_text(json.dumps(list(papers.values()),indent=2)+'\n')
 print('unique',len(rows),'priority',sum(r['screen'].startswith('PRIORITIZE') for r in rows),flush=True)

def extras(a):
 rows=json.loads((AREA/'novelty/extra_acquisition_plan.json').read_text())
 for row in rows:
  a.fetch(row['id'],row['url'],suffix=row.get('suffix','.html'),purpose=row['purpose'],max_bytes=row.get('max_bytes',LIMIT))
  print(row['id'],a.receipts[-1]['status'],flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['search','extras']);p.add_argument('--retry-unavailable',action='store_true');args=p.parse_args()
 a=Acquisition(args.mode,args.retry_unavailable);(search if args.mode=='search' else extras)(a)

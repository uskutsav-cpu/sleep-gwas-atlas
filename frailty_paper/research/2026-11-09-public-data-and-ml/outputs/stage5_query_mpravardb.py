#!/usr/bin/env python3
"""Query UCSC MPRAVarDB windows and match exact hg38 REF/ALT queue variants."""
import csv,requests,json,re,hashlib,time,zipfile,io
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
root=Path('outputs')
with (root/'stage5_score_queue.csv').open(newline='') as f:queue=list(csv.DictReader(f))
positions=defaultdict(list)
for x in queue:positions[x['chr']].append(int(x['position_GRCh38'])-1)
windows=[]
for chrom,ps in positions.items():
 start=end=None
 for pos in sorted(set(ps)):
  if start is None:start=pos;end=pos+1
  elif pos-end<=50000 and pos-start<800000:end=pos+1
  else:windows.append((chrom,start,end));start=pos;end=pos+1
 windows.append((chrom,start,end))
url='https://api.genome.ucsc.edu/getData/track'
def query(w):
 chrom,start,end=w;params={'genome':'hg38','track':'mpraVarDb','chrom':chrom,'start':start,'end':end,'maxItemsOutput':100000}
 for attempt in range(4):
  try:
   r=requests.get(url,params=params,timeout=60);r.raise_for_status()
   data=json.loads(re.sub(r'(?<=:)\s*nan(?=\s*[,}])',' null',r.text))
   if data.get('genome')!='hg38' or data.get('track')!='mpraVarDb' or 'mpraVarDb' not in data:raise ValueError('unexpected API response')
   if data.get('itemsReturned',0)>=100000 or len(data['mpraVarDb'])>=100000:raise ValueError('potential truncated response')
   return w,r.url,r.content,data['mpraVarDb']
  except Exception:
   if attempt==3:raise
   time.sleep(2**attempt)
with ThreadPoolExecutor(max_workers=3) as pool:responses=list(pool.map(query,windows))
index=defaultdict(list)
for w,u,raw,records in responses:
 for item in records:index[(item['chrom'],int(item['chromStart'])+1)].append(item)
fields=['variant_key','input_classes','input_labels','match_status','mpra_records_at_position','exact_allele_records']
matchfields=['variant_key','chrom','chromStart','chromEnd','name','ref','alt','rsid','disease','cellLine','description','log2FC','pvalue','fdr','mpraStudy','allele_orientation_note']
matches=[];coverage=[]
for q in queue:
 p=index[(q['chr'],int(q['position_GRCh38']))];exact=[v for v in p if str(v.get('ref','')).upper()==q['reference_allele'] and str(v.get('alt','')).upper()==q['alternate_allele']]
 for v in exact:matches.append({'variant_key':q['variant_key'],**{f:v.get(f,'') for f in matchfields if f not in ('variant_key','allele_orientation_note')},'allele_orientation_note':'UCSC labels match queue; verify source experiment strand, build conversion and log2FC direction before comparison'})
 coverage.append(dict(variant_key=q['variant_key'],input_classes=q['input_classes'],input_labels=q['input_labels'],match_status='exact UCSC allele labels' if exact else ('position match, allele differs' if p else 'no record in queried track'),mpra_records_at_position=len(p),exact_allele_records=len(exact)))
for name,rows,fs in [('stage5_mpravardb_coverage.csv',coverage,fields),('stage5_mpravardb_exact_allele_records.csv',matches,matchfields)]:
 with (root/name).open('w',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=fs);wr.writeheader();wr.writerows(rows)
manifest=[]
with zipfile.ZipFile(root/'stage5_mpravardb_query_responses.zip','w',zipfile.ZIP_DEFLATED) as z:
 for i,(w,url,raw,records) in enumerate(responses):
  name=f'{i:03d}_{w[0]}_{w[1]}_{w[2]}.json';z.writestr(name,raw);manifest.append({'file':name,'url':url,'sha256':hashlib.sha256(raw).hexdigest(),'records':len(records)})
(root/'stage5_mpravardb_query_receipt.json').write_text(json.dumps({'source':'UCSC MPRAVarDB hg38 track','note':'Raw API responses may contain nonstandard lowercase nan; normalized to null for parsing','queue_rows':len(queue),'windows':len(windows),'all_requests_succeeded':True,'responses':manifest},indent=2)+'\n')
print('windows',len(windows),'variants',len(queue),'exact variants',len({x['variant_key'] for x in matches}),'exact records',len(matches),'position-only',sum(x['match_status']=='position match, allele differs' for x in coverage))

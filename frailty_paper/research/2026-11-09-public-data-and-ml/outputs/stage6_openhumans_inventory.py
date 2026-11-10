#!/usr/bin/env python3
"""Aggregate public source overlap; never persist identities, URLs, or health records."""
import requests,concurrent.futures,re,json,time
from collections import Counter,defaultdict
from pathlib import Path
base='https://www.openhumans.org/api/public/datafiles/';s=requests.Session();first=s.get(base,params={'limit':1},timeout=30);first.raise_for_status();count=first.json()['count'];step=1000;offsets=list(range(0,count,step))
def fetch(offset):
 for attempt in range(3):
  try:
   r=requests.get(base,params={'limit':step,'offset':offset},timeout=40);r.raise_for_status();j=r.json();return offset,j['count'],j['results']
  except requests.RequestException:
   if attempt==2:raise
   time.sleep(1+attempt)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex: pages=list(ex.map(fetch,offsets))
seen=set();people=defaultdict(set);files=Counter();duplicates=0
for off,total,rows in pages:
 for x in rows:
  if x['id'] in seen:duplicates+=1;continue
  seen.add(x['id']);m=re.search(r'/project/(\d+)/',x.get('source_project') or '');u=x.get('user') or {}
  if not m or u.get('id') is None:continue
  src=int(m.group(1));files[src]+=1;people[src].add(u['id'])
names={128:'23andMe upload',129:'AncestryDNA upload',131:'Genome/exome upload',139:'Harvard Personal Genome Project',102:'Fitbit connection',191:'Fitbit intraday',184:'Oura Connect',453:'Apple Health port',55:'openSNP'}
selected={str(i):{'source':name,'public_files':files[i],'distinct_public_file_owners':len(people[i])} for i,name in names.items()}
genetic=people[128]|people[129]|people[131]|people[139];wearable=people[102]|people[191]|people[184]|people[453];result={'api':base,'reported_count_initial':count,'page_counts':sorted(set(total for _,total,_ in pages)),'pages':len(pages),'unique_records':len(seen),'duplicate_records_seen':duplicates,'selected_sources':selected,'genetic_source_owner_union':len(genetic),'wearable_source_owner_union':len(wearable),'intersection_owner_count':len(genetic&wearable),'caveat':'Owner overlap in live public file metadata only; source labels do not verify usable genotypes, sleep metrics, frailty outcomes, age, sample independence or complete ascertainment. Public sharing may be changed or withdrawn. PGP and other sources can contain overlapping participants.'};Path('outputs/stage6_openhumans_inventory.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

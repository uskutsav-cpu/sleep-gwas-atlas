#!/usr/bin/env python3
"""Aggregate the official project-specific public-file metadata endpoints."""
import requests,concurrent.futures,json,re,datetime
from collections import defaultdict
from pathlib import Path
names={128:'23andMe upload',129:'AncestryDNA upload',131:'Genome/exome upload',139:'Harvard Personal Genome Project',102:'Fitbit connection',191:'Fitbit intraday',184:'Oura Connect',453:'Apple Health port'}
def fetch(project):
 base=f'https://www.openhumans.org/api/public/project/{project}/datafiles/';first=requests.get(base,params={'limit':1000,'offset':0},timeout=60);first.raise_for_status();j=first.json();rows=j['results'];count=j['count'];off=len(rows);counts=[count]
 while off<count:
  r=requests.get(base,params={'limit':1000,'offset':off},timeout=60);r.raise_for_status();v=r.json();counts.append(v['count']);batch=v['results'];assert batch;rows+=batch;off+=len(batch)
 ids=set();owners=set();sleep=set();duplicates=0;wrong=0
 for x in rows:
  if x['id'] in ids:duplicates+=1;continue
  ids.add(x['id']);src=x.get('source_project','');
  if not src.rstrip('/').endswith('/'+str(project)):wrong+=1;continue
  u=x.get('user') or {};key=u.get('id');
  if key is None:continue
  owners.add(key); meta=x.get('metadata') or {};text=' '.join([str(x.get('basename','')),str(meta.get('description','')),str(meta.get('tags',''))]);
  if re.search(r'\bsleep\b',text,re.I):sleep.add(key)
 return project,{'api_reported_count':count,'page_reported_counts':sorted(set(counts)),'received_rows':len(rows),'unique_file_ids':len(ids),'duplicate_ids':duplicates,'unexpected_source':wrong,'distinct_owners':len(owners),'owners_with_sleep_keyword_metadata':len(sleep)},owners,sleep
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex: data=list(ex.map(fetch,names))
meta={str(i):{'source':names[i],**d} for i,d,_,_ in data};owners={i:s for i,_,s,_ in data};sleep={i:s for i,_,_,s in data};core=owners[128]|owners[129]|owners[131];pgp=owners[139];wear=owners[102]|owners[191]|owners[184]|owners[453];sleepwear=sleep[102]|sleep[191]|sleep[184]|sleep[453];res={'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'project_endpoints':meta,'distinct_owners_genetic_upload_projects':len(core),'distinct_owners_pgp_source':len(pgp),'distinct_owners_wearable_sources':len(wear),'distinct_owners_genetic_and_wearable':len(core&wear),'distinct_owners_pgp_and_wearable':len(pgp&wear),'distinct_owners_genetic_and_wearable_sleep_keyword':len(core&sleepwear),'limitations':'Project file metadata counts, not evidence of usable genotype, measured insomnia or frailty. Metadata keyword scan is not a data-content search; public sharing may change. People and projects may overlap.'};Path('outputs/stage6_openhumans_project_inventory.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))

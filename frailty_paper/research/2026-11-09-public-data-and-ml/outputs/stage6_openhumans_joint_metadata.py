#!/usr/bin/env python3
"""Aggregate linked metadata on explicit project endpoints; no participant IDs saved."""
from pathlib import Path
import requests,concurrent.futures,datetime,re,json
ids=[128,129,131,139,102,191,184,453]
def get(i):
 url=f'https://www.openhumans.org/api/public/project/{i}/datafiles/';r=requests.get(url,params={'limit':1000},timeout=50);r.raise_for_status();j=r.json();rows=j['results'];assert len(rows)==j['count'] and len({x['id'] for x in rows})==len(rows);return i,rows
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex: chunks=dict(ex.map(get,ids))
people={i:set() for i in ids};phenotype=set();frailty=set();sleep=set()
for i,rows in chunks.items():
 for x in rows:
  assert x['source_project'].rstrip('/').endswith('/'+str(i));k=(x.get('user') or {}).get('id');
  if k is None:continue
  people[i].add(k);meta=x.get('metadata') or {};desc=' '.join([str(x.get('basename','')),str(meta.get('description','')),str(meta.get('tags',''))])
  if i==139 and re.search(r'\b(survey|trait|phenotype|health|medical|condition)\b',desc,re.I):phenotype.add(k)
  if i==139 and re.search(r'\b(frail|frailty|insomnia)\b',desc,re.I):frailty.add(k)
  if i in (102,191,184,453) and re.search(r'\bsleep\b',desc,re.I):sleep.add(k)
g=people[128]|people[129]|people[131];w=people[102]|people[191]|people[184]|people[453];joint=g&w
result={'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'project_ids':ids,'joint_genetic_wearable_owners':len(joint),'joint_with_any_pgp_source_file':len(joint&people[139]),'joint_with_pgp_trait_keyword_metadata':len(joint&phenotype),'joint_with_pgp_frailty_or_insomnia_keyword_metadata':len(joint&frailty),'joint_with_wearable_sleep_keyword_metadata':len(joint&sleep),'interpretation':'Substring presence only in file names and tags/descriptions. Not a phenotype assessment, not proof of data content, genotype coverage, consent validity or statistical power. Absence of keyword does not imply absence of relevant measures.'}
Path('outputs/stage6_openhumans_joint_metadata.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

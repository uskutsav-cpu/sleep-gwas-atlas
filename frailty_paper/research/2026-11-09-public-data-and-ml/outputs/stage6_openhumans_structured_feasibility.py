#!/usr/bin/env python3
"""Inspect public, structured PGP trait files; emit aggregate counts only."""
from pathlib import Path
import requests,concurrent.futures,re,json,datetime,collections
ids=[128,129,131,139,102,191,184,453]
def listfiles(i):
 r=requests.get(f'https://www.openhumans.org/api/public/project/{i}/datafiles/',params={'limit':1000},timeout=50);r.raise_for_status();j=r.json();assert len(j['results'])==j['count'];return i,j['results']
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:d=dict(ex.map(listfiles,ids))
owners={i:{x['user']['id'] for x in rows if x.get('user') and x['user'].get('id')} for i,rows in d.items()};joint=(owners[128]|owners[129]|owners[131])&(owners[102]|owners[191]|owners[184]|owners[453]);selected=[]
for x in d[139]:
 u=(x.get('user') or {}).get('id');m=x.get('metadata') or {};name=x.get('basename') or '';desc=' '.join([name,str(m.get('description','')),str(m.get('tags',''))]);
 if u in joint and name.lower().endswith('.json') and re.search(r'\b(survey|trait|phenotype|health|medical|condition)\b',desc,re.I):selected.append(x)
patterns={'insomnia':r'\binsomnia\b','frailty':r'\bfrail(?:ty)?\b','sleep':r'\bsleep\b','fatigue':r'\b(fatigue|exhaustion|exhausted)\b','walking':r'\b(walk|walking|gait)\b','grip':r'\b(grip|handgrip)\b','weight_loss':r'\bweight\s+loss\b','physical_activity':r'\bphysical\s+activity\b'}
found={k:set() for k in patterns};ok=set();errors=collections.Counter();n_valid=0
for x in selected:
 u=x['user']['id'];url=x.get('download_url','')
 try:
  if not url.startswith('https://'):errors['non_https']+=1;continue
  with requests.get(url,stream=True,timeout=25) as r:
   if not r.ok:errors['http_failure']+=1;continue
   buf=bytearray()
   for chunk in r.iter_content(65536):
    buf.extend(chunk)
    if len(buf)>5_000_000:break
   if len(buf)>5_000_000:errors['size_limit']+=1;continue
  obj=json.loads(buf);text=json.dumps(obj,ensure_ascii=False);n_valid+=1;ok.add(u)
  for key,pattern in patterns.items():
   if re.search(pattern,text,re.I):found[key].add(u)
 except (requests.RequestException,ValueError,UnicodeError):errors['request_or_json_error']+=1
out={'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'joint_owners':len(joint),'selected_json_files':len(selected),'valid_json_files':n_valid,'owners_with_valid_json':len(ok),'errors':dict(errors),'owners_with_term_anywhere':{k:len(v) for k,v in found.items()},'limitations':'Text keyword in unvalidated public JSON, not a patient diagnosis or a constructed phenotype; no values or identities retained. No inspection of non-JSON records, genotype files, or outcomes. Public sharing can change.'};Path('outputs/stage6_openhumans_structured_feasibility.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

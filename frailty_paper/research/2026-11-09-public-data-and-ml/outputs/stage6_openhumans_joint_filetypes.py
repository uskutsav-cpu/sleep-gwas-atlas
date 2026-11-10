from pathlib import Path
import requests,concurrent.futures,re,json,collections,datetime
ids=[128,129,131,139,102,191,184,453]
def get(i):
 r=requests.get(f'https://www.openhumans.org/api/public/project/{i}/datafiles/',params={'limit':1000},timeout=50);r.raise_for_status();j=r.json();assert len(j['results'])==j['count'];return i,j['results']
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:d=dict(ex.map(get,ids))
owners={i:{x['user']['id'] for x in rows if x.get('user') and x['user'].get('id')} for i,rows in d.items()};joint=(owners[128]|owners[129]|owners[131])&(owners[102]|owners[191]|owners[184]|owners[453]);types=collections.Counter();owners_types=collections.defaultdict(set);total=0
for x in d[139]:
 user=(x.get('user') or {}).get('id');
 if user not in joint:continue
 meta=x.get('metadata') or {};t=' '.join([str(x.get('basename','')),str(meta.get('description','')),str(meta.get('tags',''))]);
 if not re.search(r'\b(survey|trait|phenotype|health|medical|condition)\b',t,re.I):continue
 total+=1;name=(x.get('basename') or '').lower();kind='json' if name.endswith('.json') else 'csv' if name.endswith('.csv') else 'zip' if name.endswith('.zip') else 'pdf' if name.endswith('.pdf') else 'text' if name.endswith(('.txt','.tsv')) else 'other';types[kind]+=1;owners_types[kind].add(user)
r={'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'matched_owners':len(joint),'pgp_trait_keyword_files':total,'file_types':{k:{'files':v,'distinct_owners':len(owners_types[k])} for k,v in types.items()},'caveat':'Metadata-only filename extensions and keyword matches. No assessment of file content, current download access, clinical phenotype or genetic variant coverage; owners can occur in multiple types.'};Path('outputs/stage6_openhumans_joint_filetypes.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

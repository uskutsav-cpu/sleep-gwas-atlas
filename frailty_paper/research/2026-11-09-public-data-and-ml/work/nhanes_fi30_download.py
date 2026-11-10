import requests,pathlib,concurrent.futures,hashlib,json
jobs=[(y,s,n) for y,s in [(2011,'G'),(2013,'H')] for n in ['CDQ','DPQ','BPQ','MCQ','ENQ','DIQ','OSQ','PFD','OHQ','PFQ']]
def get(j):
 y,s,n=j;u=f'https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{y}/DataFiles/{n}_{s}.XPT';r=requests.get(u,timeout=60);p=pathlib.Path(f'work/nhanes{y}/{n}_{s}.XPT');d={'year':y,'name':n,'status':r.status_code,'url':u}
 if r.ok and r.content[:6]==b'HEADER':p.write_bytes(r.content);d.update(bytes=len(r.content),sha256=hashlib.sha256(r.content).hexdigest())
 return d
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as e:out=list(e.map(get,jobs))
open('outputs/stage6_nhanes_fi30_sources.json','w').write(json.dumps(out,indent=2)+'\n');print([(x['year'],x['name'],x['status'],x.get('bytes')) for x in out])

import hashlib,json
from pathlib import Path
import requests
p=Path(__file__).resolve().parent
urls={'record':'https://zenodo.org/api/records/14011550','readme':'https://zenodo.org/records/14011550/files/0_Readme.txt?download=1','repository_tree':'https://api.github.com/repos/JiajinChen/PRS-HF/git/trees/HEAD?recursive=1'}
objects={}
for key,u in urls.items():
 r=requests.get(u,timeout=60);r.raise_for_status();(p/(key+('.txt' if key=='readme' else '.json'))).write_bytes(r.content);objects[key]={'url':u,'sha256':hashlib.sha256(r.content).hexdigest(),'bytes':len(r.content)}
z=json.loads((p/'record.json').read_text());t=json.loads((p/'repository_tree.json').read_text());read=(p/'readme.txt').read_text();assert z['metadata']['creators'][0]['name']=='Chen, Jiajin' and '47139' in read and '349680' in read
objects['repository_tree'].update(commit=t['sha'],truncated=t['truncated'],paths=len(t['tree']));(p/'receipt.json').write_text(json.dumps(objects,indent=2)+'\n');print(json.dumps(objects,indent=2))

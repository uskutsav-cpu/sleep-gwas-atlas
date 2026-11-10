from pathlib import Path
import requests,subprocess,json,hashlib,re
O=Path('outputs');W=Path('work/pgp_survey_forms');W.mkdir(exist_ok=True)
url='https://api.github.com/repos/PGPHarvard/pgp-surveys/contents/Surveys';r=requests.get(url,timeout=30);r.raise_for_status();files=[x for x in r.json() if x['name'].endswith('.pdf')];results=[]
for i,x in enumerate(files):
 b=requests.get(x['download_url'],timeout=30);b.raise_for_status();z=subprocess.run(['pdftotext','-layout','-','-'],input=b.content,capture_output=True,check=True).stdout.decode('utf-8','replace'); matches=[]
 for line in z.splitlines():
  if re.search(r'\b(insomnia|sleep|frail|walk|walking|grip|weakness|fatigue|exhaust|weight loss|activity)\b',line,re.I):matches.append(line.strip())
 results.append({'file':x['name'],'url':x['html_url'],'sha256':hashlib.sha256(b.content).hexdigest(),'extracted_characters':len(z),'matches':matches})
(O/'stage6_pgp_survey_inventory.json').write_text(json.dumps({'source':url,'files':results},indent=2)+'\n');print('forms',len(results))
for x in results:
 if x['matches']: print(x['file'],':',' | '.join(x['matches'][:12]))

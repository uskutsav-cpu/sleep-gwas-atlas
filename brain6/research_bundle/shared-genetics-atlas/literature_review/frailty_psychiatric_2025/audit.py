import hashlib,json,subprocess,requests
from pathlib import Path
p=Path(__file__).resolve().parent
r=requests.get('https://zenodo.org/api/records/15449098',timeout=30);r.raise_for_status();z=r.json();item=z['files'][0];sources={'zenodo':item['links']['self'],'university':'https://www.symc.edu.cn/lw/1guohongGenome.pdf'};out={}
for k,u in sources.items():

 try:r=requests.get(u,timeout=90);r.raise_for_status()
 except requests.exceptions.RequestException as e:
  out[k]={'url':u,'status':'download_failed','error_type':type(e).__name__};continue
 data=r.content;assert data.startswith(b'%PDF-');proc=subprocess.run(['pdftotext','-layout','-','-'],input=data,capture_output=True,check=True);text=proc.stdout.decode('utf-8',errors='replace');out[k]={'url':u,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'text_sha256':hashlib.sha256(proc.stdout).hexdigest(),'page_breaks':text.count('\f'),'title_present':'Genome-Wide Crossover Analysis' in text,'frailty_source_present':'GCST90020053' in text,'finn_sleep_present':'finngen_R11_SLEEP' in text,'journal_header_machine_learning':'Machine Learning' in text,'journal_header_organizational':'International Journal of Interdisciplinary Organizational Studies' in text};print(k,len(data),out[k]['sha256'])
q='Genome-Wide Crossover Analysis Identification of Shared Genetic Effect Site between Frailty and Multiple Psychiatric Disorders';r=requests.get('https://api.crossref.org/works',params={'query.title':q,'rows':10},timeout=30);r.raise_for_status();items=r.json()['message']['items'];exact=[{'doi':x.get('DOI'),'title':x.get('title',[None])[0]} for x in items if x.get('title') and x['title'][0].strip().casefold()==q.casefold()];info={'zenodo_doi':z['metadata'].get('doi'),'zenodo_deposit_type':z['metadata'].get('resource_type'),'zenodo_date':z['metadata'].get('publication_date'),'files':out,'crossref_query_title':q,'crossref_exact_title_matches_top_ten':exact,'note':'Search is limited; absence in these ten returned matches is not proof of nonpublication. PDF header identity requires independent verification.'};(p/'receipt.json').write_text(json.dumps(info,indent=2)+'\n');print('exact Crossref matches:',len(exact))

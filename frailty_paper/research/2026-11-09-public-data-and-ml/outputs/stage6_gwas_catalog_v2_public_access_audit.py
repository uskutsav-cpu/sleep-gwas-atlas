import requests,json,concurrent.futures,hashlib
from pathlib import Path
u='https://www.ebi.ac.uk/gwas/rest/api/v2/studies';pmids=['35670608','37523034','41733842','42020437','37184129'];out=Path('outputs')
def f(pmid):
 r=requests.get(u,params={'pubmed_id':pmid,'size':100},timeout=45);r.raise_for_status();j=r.json();page=j.get('page',{});assert page.get('totalPages',1)<=1;body=r.content;(out/f'stage6_gwas_catalog_v2_pmid_{pmid}.json').write_bytes(body);studies=j.get('_embedded',{}).get('studies',[]);return {'pmid':pmid,'url':r.url,'total':page.get('totalElements'),'sha256':hashlib.sha256(body).hexdigest(),'studies':[{'accession':x.get('accession_id'),'trait':x.get('disease_trait'),'full_summary_stats_available':x.get('full_summary_stats_available'),'full_summary_stats':x.get('full_summary_stats')} for x in studies]}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:x=list(ex.map(f,pmids))
(out/'stage6_gwas_catalog_v2_public_access_audit.json').write_text(json.dumps(x,indent=2)+'\n');print(json.dumps(x,indent=2))

#!/usr/bin/env python3
"""Bounded public metadata acquisition; no GWAS downloads or access bypass."""
import hashlib,json,re,time,urllib.parse,zipfile,subprocess
from pathlib import Path
from xml.etree import ElementTree as ET
from audit_evidence import ROOT,AREA,write_tsv

STORE=ROOT/'work/literature-20261007'
MAX_BYTES=30_000_000
HEADERS={'User-Agent':'Brain6ScientificAudit/1.0 (public literature metadata only)'}
QUERY_NAMES={
 'sleep_psychiatric':'(sleep OR insomnia OR chronotype) AND (psychiatric OR depression OR schizophrenia OR bipolar OR ADHD) AND (GWAS OR genetic OR colocalization OR pleiotropic)',
 'sleep_neurodegenerative':'(sleep OR insomnia OR chronotype) AND (Parkinson OR Alzheimer OR neurodegenerative) AND (GWAS OR genetic OR colocalization)',
 'local_methods':'(sleep OR insomnia) AND (LAVA OR SUPERGNOVA OR "local genetic" OR "fine mapping")',
 'sleep_gwas':'(sleep OR insomnia OR chronotype) AND ("genome-wide association" OR GWAS)',
 'preprints':'(sleep OR insomnia) AND (genetic OR GWAS OR colocalization) AND (SRC:PPR)'}

def main():
 STORE.mkdir(parents=True,exist_ok=True); manifest=[];receipts=[];papers={};tables=[]
 for name,q in QUERY_NAMES.items():
  q+=' AND FIRST_PDATE:[2025-01-01 TO 2026-10-07]'
  url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode({'query':q,'format':'json','resultType':'core','pageSize':1000})
  manifest.append({'id':name,'url':url,'max_bytes':MAX_BYTES,'purpose':'Date-bounded literature search; no absence-based novelty claim','access':'PUBLIC_METADATA','query':q})
 for name,url in [
  ('jia_article','https://academic.oup.com/sleep/article/48/1/zsae209/7750695'),
  ('xue_article','https://academic.oup.com/sleep/article/49/1/zsaf317/8279895'),
  ('zu_article','https://www.nature.com/articles/s41398-026-04166-4'),
  ('genome_guidelines','https://link.springer.com/journal/13073/submission-guidelines/research'),
  ('genome_collections','https://link.springer.com/journal/13073/collections'),
  ('mp_guidelines','https://www.nature.com/mp/authors-and-referees/preparation-of-articles'),
  ('gwas_api_docs','https://www.ebi.ac.uk/gwas/docs/programmatic-access/rest-api/'),
  ('finngen_access','https://www.finngen.fi/en/access_results')]:
  manifest.append({'id':name,'url':url,'max_bytes':MAX_BYTES,'purpose':'Inspect source text and actual supplement links','access':'PUBLIC_WEB','query':''})
 write_tsv(AREA/'review/acquisition_manifest.tsv',manifest)
 def fetch(item):
  record=dict(item);record['retrieved_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
  try:
   temporary=STORE/(item['id']+'.download')
   result=subprocess.run(['curl','--fail','--location','--silent','--show-error','--max-time','30',
      '--max-filesize',str(MAX_BYTES),'--user-agent',HEADERS['User-Agent'],'--output',str(temporary),
      '--write-out','%{json}',item['url']],capture_output=True,text=True)
   if result.returncode:raise RuntimeError(result.stderr.strip())
   metadata=json.loads(result.stdout);payload=temporary.read_bytes()
   record.update(http_status=metadata['http_code'],final_url=metadata['url_effective'],content_type=metadata.get('content_type',''))
   if len(payload)>MAX_BYTES:raise ValueError('object exceeds bounded acquisition manifest')
   suffix='.json' if 'json' in record['content_type'] else '.docx' if payload[:2]==b'PK' else '.html'
   path=STORE/(item['id']+suffix);temporary.rename(path)
   record.update(status='RETRIEVED',bytes=len(payload),sha256=hashlib.sha256(payload).hexdigest(),local_path=str(path.relative_to(ROOT)))
   receipts.append(record);return payload,path
  except Exception as e:
   record.update(status='UNAVAILABLE',error=str(e));receipts.append(record);return None,None
 def checkpoint():
  (AREA/'review/acquisition_receipts.json').write_text(json.dumps(receipts,indent=2,sort_keys=True)+'\n')
 for item in manifest:
  payload,path=fetch(item);checkpoint()
  if not payload:continue
  if item['id'] in QUERY_NAMES:
   d=json.loads(payload);print(item['id'],'hits',d['hitCount'],'retrieved',len(d['resultList']['result']))
   for r in d['resultList']['result']:
    key=(r.get('source',''),r.get('id',''));papers[key]=r
   if d['hitCount']>len(d['resultList']['result']):
    print('TRUNCATED_SEARCH_REQUIRES_PAGINATION',item['id'])
  elif item['id'] in {'jia_article','xue_article','zu_article'}:
   html=payload.decode('utf-8',errors='replace')
   links=re.findall(r'(?:href|data-track-action)=["\']([^"\']+)["\']',html)
   supplements=[]
   for link in links:
    link=link.replace('&amp;','&')
    if re.search(r'\.(docx|xlsx|zip)(\?|$)',link,re.I) and ('supp' in link.lower() or 'MOESM' in link):
     url=urllib.parse.urljoin(receipts[-1]['final_url'],link)
     if url not in supplements:supplements.append(url)
   for i,url in enumerate(supplements):
    si={'id':item['id']+f'_supp{i+1}','url':url,'max_bytes':MAX_BYTES,'purpose':'Source-linked supplementary locus inspection','access':'PUBLIC_LINK_FROM_RETRIEVED_ARTICLE','query':''}
    manifest.append(si);write_tsv(AREA/'review/acquisition_manifest.tsv',manifest)
    data,sp=fetch(si);checkpoint()
    if data and data[:2]==b'PK' and zipfile.is_zipfile(sp):
     with zipfile.ZipFile(sp) as z:
      if 'word/document.xml' in z.namelist():
       doc=ET.fromstring(z.read('word/document.xml'));ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
       paras=[''.join(n.itertext()) for n in doc.findall('.//w:p',ns)]
       (STORE/(si['id']+'.txt')).write_text('\n'.join(paras))
       for ti,table in enumerate(doc.findall('.//w:tbl',ns)):
        rr=[]
        for tr in table.findall('./w:tr',ns):
         rr.append([''.join(t.text or '' for t in cell.findall('.//w:t',ns)) for cell in tr.findall('./w:tc',ns)])
        tables.append({'source':si['id'],'table_index':ti+1,'rows':rr})
 for key,r in sorted(papers.items()):
  # Archive metadata privately; derived bibliography and selection log are committed.
  (STORE/('paper_'+key[0]+'_'+key[1]+'.json')).write_text(json.dumps(r,indent=2)+'\n')
 selected=[{'source':r.get('source',''),'id':r.get('id',''),'pmid':r.get('pmid',''),'doi':r.get('doi',''),'first_publication_date':r.get('firstPublicationDate',''),'title':r.get('title',''),'authors':r.get('authorString',''),'journal':r.get('journalTitle',''),'pmcid':r.get('pmcid',''),'review_status':'TITLE_METADATA_RETRIEVED_NOT_FULLTEXT_REVIEWED'} for r in papers.values()]
 write_tsv(AREA/'review/literature_search_records.tsv',selected,fields=['source','id','pmid','doi','first_publication_date','title','authors','journal','pmcid','review_status'])
 (STORE/'supplement_tables.json').write_text(json.dumps(tables,indent=2)+'\n')
 checkpoint();print('unique records',len(papers),'supplement tables',len(tables))

if __name__=='__main__':main()

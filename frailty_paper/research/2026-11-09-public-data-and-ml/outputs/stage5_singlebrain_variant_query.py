#!/usr/bin/env python3
"""Retrieve public SingleBrain portal summaries for nominated genes and match exact GRCh38 alleles."""
from pathlib import Path
import concurrent.futures, csv, json, re, time, zipfile
import requests
BASE='https://singlebrain.nygenome.org/eqtl/fetch_sum.php'
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'outputs'
genes={'SMAD5','FOXP2','INKA1','TMEM161B-DT'}
top=list(csv.DictReader(open(OUT/'stage5_singlebrain_gene_top_associations.csv',newline='')))
queue={r['variant_key']:r for r in csv.DictReader(open(OUT/'stage5_score_queue.csv',newline=''))}
tasks=sorted({(r['cell_type'],r['gene_symbol'],r['feature']) for r in top if r['gene_symbol'] in genes})
def query(task):
 cell,gene,feature=task;p={'phenotype_id':feature,'condition_id':cell,'draw':'1','start':'0','length':'10000','search[value]':'','search[regex]':'false','order[0][column]':'4','order[0][dir]':'desc'}
 for i in range(5):
  for k,v in {'data':str(i),'name':'','searchable':'true','orderable':'true'}.items():p[f'columns[{i}][{k}]']=v
  p[f'columns[{i}][search][value]']='';p[f'columns[{i}][search][regex]']='false'
 for attempt in range(3):
  try:
   r=requests.get(BASE,params=p,timeout=45);r.raise_for_status();j=r.json()
   if j['recordsFiltered']!=len(j['data']):raise ValueError('Incomplete response: '+str((j['recordsFiltered'],len(j['data']))))
   return task,j
  except Exception:
   if attempt==2:raise
   time.sleep(2*(attempt+1))
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool: results=list(pool.map(query,tasks))
rows=[]; receipt=[]
with zipfile.ZipFile(OUT/'stage5_singlebrain_priority_gene_portal_raw.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
 for (cell,gene,feature),j in results:
  z.writestr(f'{cell}_{gene}_{feature}.json',json.dumps(j,separators=(',',':')))
  matched=0
  for rsid,posid,beta,p,zscore in j['data']:
   m=re.fullmatch(r'(?:chr)?(\d+|X|Y):(\d+)_([ACGT]+)\/([ACGT]+)',str(posid),re.I)
   if not m:continue
   chrom,pos,ref,alt=m.groups(); key=f'chr{chrom}:{pos}:{ref.upper()}>{alt.upper()}'
   if key not in queue:continue
   matched+=1; source=queue[key]
   rows.append({'cell_type':cell,'gene':gene,'feature':feature,'variant_key':key,'rsid':rsid,'portal_beta':beta,'portal_random_p':p,'portal_z_abs_or_signed':zscore,'input_classes':source['input_classes'],'gene_context':source['gene_context']})
  receipt.append({'cell_type':cell,'gene':gene,'feature':feature,'records_total':j['recordsTotal'],'records_returned':len(j['data']),'exact_queue_matches':matched})
fields=['cell_type','gene','feature','variant_key','rsid','portal_beta','portal_random_p','portal_z_abs_or_signed','input_classes','gene_context']
with open(OUT/'stage5_singlebrain_priority_variant_matches.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
(OUT/'stage5_singlebrain_priority_query_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('gene-cell queries',len(results),'total records',sum(i['records_returned'] for i in receipt),'exact match rows',len(rows),'distinct variants',len({r['variant_key'] for r in rows}))
for i in receipt:print(i['cell_type'],i['gene'],i['records_returned'],i['exact_queue_matches'])

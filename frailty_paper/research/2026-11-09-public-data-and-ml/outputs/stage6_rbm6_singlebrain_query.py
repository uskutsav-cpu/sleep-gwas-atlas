#!/usr/bin/env python3
from pathlib import Path
import concurrent.futures,requests,json,zipfile,csv,re
root=Path(__file__).resolve().parents[1];out=root/'outputs'; j=json.loads((out/'stage6_rbm6_five_gwas_region_summary.json').read_text());variants={v['rsid']:v for x in j['comparisons'].values() for v in [x['sleep_best_shared_variant'],x['frailty_best_shared_variant']]};variants['rs2353586']={'rsid':'rs2353586','pos':50008613,'ref':'G','alt':'T'}
cells=['Ast','End','Ext','IN','MG','OD','OPC']
def query(cell):
 p={'phenotype_id':'ENSG00000004534.15','condition_id':cell,'draw':'1','start':'0','length':'10000','search[value]':'','search[regex]':'false','order[0][column]':'4','order[0][dir]':'desc'}
 for i in range(5):
  for k,v in {'data':str(i),'name':'','searchable':'true','orderable':'true'}.items():p[f'columns[{i}][{k}]']=v
  p[f'columns[{i}][search][value]']='';p[f'columns[{i}][search][regex]']='false'
 r=requests.get('https://singlebrain.nygenome.org/eqtl/fetch_sum.php',params=p,timeout=45);r.raise_for_status();data=r.json();assert len(data['data'])==data['recordsFiltered'];return cell,data
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:results=list(ex.map(query,cells))
rows=[];coverage=[]
with zipfile.ZipFile(out/'stage6_rbm6_singlebrain_raw.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
 for cell,data in results:
  z.writestr(cell+'.json',json.dumps(data,separators=(',',':')));count=0
  for rsid,posid,beta,p,zscore in data['data']:
   if rsid not in variants:continue
   v=variants[rsid];expected=f"3:{v['pos']}_{v['ref']}/{v['alt']}"
   if posid!=expected:continue
   count+=1;rows.append({'cell_type':cell,'gene':'RBM6','gene_id':'ENSG00000004534.15','variant':rsid,'variant_GRCh38':f"chr3:{v['pos']}:{v['ref']}>{v['alt']}",'portal_beta':beta,'portal_random_p':p,'portal_z':zscore})
  coverage.append({'cell_type':cell,'records_returned':len(data['data']),'exact_matches':count})
with open(out/'stage6_rbm6_singlebrain_exact_matches.csv','w',newline='') as f:
 fields=['cell_type','gene','gene_id','variant','variant_GRCh38','portal_beta','portal_random_p','portal_z'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
(out/'stage6_rbm6_singlebrain_query_receipt.json').write_text(json.dumps(coverage,indent=2)+'\n');print('exact rows',len(rows));print('\n'.join(str(x) for x in rows if x['variant']=='rs2353586'))

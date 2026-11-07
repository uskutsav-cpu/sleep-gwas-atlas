import csv,json,pathlib,hashlib,subprocess,urllib.parse
p=pathlib.Path(__file__).resolve().parent;targets=[('rs817509','chr3_118032014_A_G_b38','ENSG00000185565.11'),('rs12536335','chr7_114403104_A_G_b38','ENSG00000128573.24')];out=[]
for rs,vid,gene in targets:
 for tissue in ['Brain_Cortex','Brain_Frontal_Cortex_BA9']:
  url='https://gtexportal.org/api/v2/association/dyneqtl?'+urllib.parse.urlencode({'datasetId':'gtex_v8','variantId':vid,'gencodeId':gene,'tissueSiteDetailId':tissue})
  try:
   z=subprocess.run(['curl','--fail','--silent','--show-error','--max-time','30',url],capture_output=True,check=True).stdout;d=json.loads(z);assert d['variantId']==vid and d['gencodeId']==gene and d['tissueSiteDetailId']==tissue
   n=sum(d[k] for k in ['hetCount','homoAltCount','homoRefCount']);out.append({'status':'OK','variant':rs,'variant_id':vid,'gene':d['geneSymbol'],'gencode_id':gene,'tissue':tissue,'n':n,'nominal_p':d['pValue'],'nes':d['nes'],'se':d['error'],'gene_threshold':d.get('pValueThreshold'),'url':url,'response_sha256':hashlib.sha256(z).hexdigest()})
  except Exception as e:out.append({'status':'ERROR '+str(e),'variant':rs,'variant_id':vid,'gencode_id':gene,'tissue':tissue,'url':url})
  print(rs,tissue,out[-1]['status'],out[-1].get('nominal_p'))
fields=list(dict.fromkeys(k for r in out for k in r))
with (p/'gtex_dynamic_brain_summary.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(out)

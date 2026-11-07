import csv,json,pathlib,subprocess,urllib.parse,hashlib
out=pathlib.Path(__file__).resolve().parent; variant='chr7_114693988_T_C_b38'; gene='ENSG00000128573.24'; base='https://gtexportal.org/api/v2/association/';status=[];associations=[]
def fetch(url):
 r=subprocess.run(['curl','--fail','--silent','--show-error','--max-time','45',url],check=True,capture_output=True);return r.stdout,json.loads(r.stdout)
for typ in ['Eqtl','Sqtl']:
 params={'datasetId':'gtex_v8','variantId':variant,'itemsPerPage':1000,'page':0};url=base+'singleTissue'+typ+'?'+urllib.parse.urlencode(params)
 try:
  raw,j=fetch(url); assert isinstance(j.get('data'),list); data=j['data']; total=j.get('total',len(data)); assert total==len(data)
  assert all(x.get('variantId')==variant for x in data);(out/f'{typ}_response.json').write_bytes(raw)
  status.append({'type':typ,'status':'OK','n':len(data),'url':url,'sha256':hashlib.sha256(raw).hexdigest()})
  for x in data:associations.append({'type':typ,**{k:x.get(k,'') for k in ('snpId','variantId','gencodeId','geneSymbol','tissueSiteDetailId','pValue','nes','phenotypeId')}})
 except Exception as e:status.append({'type':typ,'status':'ERROR '+str(e),'n':'','url':url,'sha256':''})
for tissue in ['Brain_Cortex','Brain_Frontal_Cortex_BA9','Adrenal_Gland']:
 params={'datasetId':'gtex_v8','tissueSiteDetailId':tissue,'gencodeId':gene,'variantId':variant};url=base+'dyneqtl?'+urllib.parse.urlencode(params)
 try:
  raw,j=fetch(url);assert j['variantId']==variant and j['gencodeId']==gene and j['tissueSiteDetailId']==tissue
  status.append({'type':'dynamic_'+tissue,'status':'OK','n':sum(j[k] for k in ['hetCount','homoAltCount','homoRefCount']),'url':url,'sha256':hashlib.sha256(raw).hexdigest(),'gene':j.get('geneSymbol'),'pValue':j.get('pValue'),'nes':j.get('nes'),'se':j.get('error'),'pValueThreshold':j.get('pValueThreshold')})
 except Exception as e:status.append({'type':'dynamic_'+tissue,'status':'ERROR '+str(e),'n':'','url':url,'sha256':''})
def write(name,rows,fields):
 with (out/name).open('w',newline='') as f:w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)
write('query_status.tsv',status,['type','status','n','url','sha256','gene','pValue','nes','se','pValueThreshold'])
write('significant_pairs.tsv',associations,['type','snpId','variantId','gencodeId','geneSymbol','tissueSiteDetailId','pValue','nes','phenotypeId'])
for r in status:print(r)

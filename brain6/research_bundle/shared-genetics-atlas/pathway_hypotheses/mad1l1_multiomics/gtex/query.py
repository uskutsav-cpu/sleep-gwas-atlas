import pathlib,subprocess,urllib.parse,json,csv,hashlib
p=pathlib.Path(__file__).resolve().parent;raw=p/'responses';raw.mkdir(exist_ok=True)
variants={'rs28528997':'chr7_2034487_A_G_b38','rs56305291':'chr7_2009071_A_G_b38'}
status=[];assoc=[]
for rs,var in variants.items():
 for typ in ('Eqtl','Sqtl'):
  url='https://gtexportal.org/api/v2/association/singleTissue'+typ+'?'+urllib.parse.urlencode({'datasetId':'gtex_v8','variantId':var,'itemsPerPage':1000,'page':0})
  item={'rsid':rs,'variant_id':var,'type':typ,'url':url}
  try:
   proc=subprocess.run(['curl','--fail','--silent','--show-error','--max-time','40',url],capture_output=True,check=True);data=proc.stdout;j=json.loads(data);arr=j['data'];total=int(j['paging_info']['totalNumberOfItems']);assert total==len(arr) and total<=1000 and all(x['variantId']==var and x['datasetId']=='gtex_v8' for x in arr)
   dest=raw/(rs+'_'+typ+'.json');dest.write_bytes(data);h=hashlib.sha256(data).hexdigest();item.update(status='OK',n=len(arr),sha256=h)
   for x in arr:assoc.append({'rsid':rs,'variant_id':var,'type':typ,'gene':x.get('geneSymbol',''),'gencode_id':x.get('gencodeId',''),'tissue':x.get('tissueSiteDetailId',''),'p':x.get('pValue',''),'nes':x.get('nes',''),'phenotype_id':x.get('phenotypeId',''),'response_sha256':h})
  except Exception as e:item.update(status='ERROR '+str(e))
  status.append(item);print(rs,typ,item['status'],item.get('n',''),flush=True)
for name,data,fields in [('queries.tsv',status,['rsid','variant_id','type','status','n','url','sha256']),('associations.tsv',assoc,['rsid','variant_id','type','gene','gencode_id','tissue','p','nes','phenotype_id','response_sha256'])]:
 with (p/name).open('w',newline='') as f:w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(data)

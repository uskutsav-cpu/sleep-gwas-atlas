import csv,json,pathlib,hashlib,subprocess,urllib.parse
p=pathlib.Path(__file__).resolve().parent;variants={'rs817509':'chr3_118032014_A_G_b38','rs12536335':'chr7_114403104_A_G_b38'};q=[];rows=[]
for rs,vid in variants.items():
 for typ in ['Eqtl','Sqtl']:
  url='https://gtexportal.org/api/v2/association/singleTissue'+typ+'?'+urllib.parse.urlencode({'datasetId':'gtex_v8','variantId':vid,'itemsPerPage':1000});z=subprocess.run(['curl','--fail','--silent','--show-error','--max-time','30',url],capture_output=True,check=True).stdout;d=json.loads(z);assert len(d['data'])==d['paging_info']['totalNumberOfItems']<=1000 and all(x['variantId']==vid for x in d['data']);f=p/'gtex_responses'/(rs+'_'+typ+'.json');f.write_bytes(z);h=hashlib.sha256(z).hexdigest();q.append({'variant':rs,'type':typ,'variant_id':vid,'n':len(d['data']),'url':url,'response_sha256':h})
  for r in d['data']:rows.append({'variant':rs,'type':typ,'variant_id':vid,'gene':r.get('geneSymbol',''),'gencode_id':r.get('gencodeId',''),'tissue':r.get('tissueSiteDetailId',''),'nominal_p':r.get('pValue',''),'nes':r.get('nes',''),'splice_phenotype_id':r.get('phenotypeId','')})
  print(rs,typ,len(d['data']))
for name,arr in [('gtex_query_receipt.tsv',q),('gtex_associations.tsv',rows)]:
 with (p/name).open('w',newline='') as f:
  cols=list(arr[0]) if arr else ['variant'];w=csv.DictWriter(f,cols,delimiter='\t');w.writeheader();w.writerows(arr)

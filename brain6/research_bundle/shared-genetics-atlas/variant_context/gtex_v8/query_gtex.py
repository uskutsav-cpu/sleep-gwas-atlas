import csv,json,hashlib,subprocess,pathlib,time,urllib.parse
root=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/brain6/outputs/shared-genetics-atlas/variant_context');out=root/'gtex_v8'; raw=out/'responses';raw.mkdir(parents=True,exist_ok=True)
with (root/'selected_variant_gene_context.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
queries=[]; assoc=[]
for r in rows:
 chrom,pos,alleles=r['grch38_ref_alt'].split(':');ref,alt=alleles.split('>'); var=f'{chrom}_{pos}_{ref}_{alt}_b38'
 for kind in ['Eqtl','Sqtl']:
  params={'datasetId':'gtex_v8','variantId':var,'itemsPerPage':1000,'page':0};url='https://gtexportal.org/api/v2/association/singleTissue'+kind+'?'+urllib.parse.urlencode(params);path=raw/(r['snp']+'_'+kind+'.json')
  try:
   if path.exists():data=path.read_bytes()
   else:
    data=subprocess.run(['curl','--fail','--silent','--show-error','--max-time','30',url],check=True,capture_output=True).stdout
    parsed=json.loads(data);path.write_bytes(data);time.sleep(.1)
   j=json.loads(data); arr=j['data'];paging=j['paging_info'];assert isinstance(arr,list);total=int(paging['totalNumberOfItems']);assert total<=1000 and total==len(arr);assert all(x['variantId']==var and x['datasetId']=='gtex_v8' for x in arr)
   queries.append({'snp':r['snp'],'variant_id':var,'type':kind,'status':'OK','n_returned':len(arr),'total_reported':total,'url':url,'sha256':hashlib.sha256(data).hexdigest(),'local_response':str(path.relative_to(out))})
   for x in arr:
    assoc.append({'snp':r['snp'],'type':kind,'variant_id':var,'gene_symbol':x.get('geneSymbol',''),'gencode_id':x.get('gencodeId',''),'tissue':x.get('tissueSiteDetailId',''),'p_value':x.get('pValue',''),'normalized_effect_size':x.get('nes',''),'phenotype_id':x.get('phenotypeId',''),'query_sha256':hashlib.sha256(data).hexdigest()})
   print(r['snp'],kind,len(arr),flush=True)
  except Exception as e:
   queries.append({'snp':r['snp'],'variant_id':var,'type':kind,'status':'ERROR: '+str(e),'url':url})
   print(r['snp'],kind,'ERROR',e,flush=True)
for name,data,fields in [('query_status.tsv',queries,['snp','variant_id','type','status','n_returned','total_reported','url','sha256','local_response']),('significant_associations.tsv',assoc,['snp','type','variant_id','gene_symbol','gencode_id','tissue','p_value','normalized_effect_size','phenotype_id','query_sha256'])]:
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(data)
print('queries',len(queries),'errors',sum(x['status']!='OK' for x in queries),'associations',len(assoc))

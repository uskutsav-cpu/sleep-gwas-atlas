import csv,json,pathlib,subprocess,concurrent.futures,hashlib,urllib.parse,time
root=pathlib.Path('.'); source=list(csv.DictReader(open('outputs/shared22_eight_source_evidence_matrix.csv'))); extra=next(x for x in csv.DictReader(open('outputs/integrated_15_candidate_review.csv')) if x['variant']=='rs323509')
targets=[dict(variant=x['variant'],chrom=x['chr_GRCh38'],pos=x['position_GRCh38'],a=x['score_counted_allele'],b=x['other_allele']) for x in source];targets.append(dict(variant='rs323509',chrom=extra['GRCh37_chr'],pos=extra['GRCh38_pos'],a=extra['counted_allele'],b=extra['other_allele']));assert len(targets)==23
raw=root/'work/stage3_gtex_v10_raw';raw.mkdir(exist_ok=True)
def query(task):
 x,kind,ref,alt=task; variantid=f"chr{x['chrom']}_{x['pos']}_{ref}_{alt}_b38"; path=f"https://gtexportal.org/api/v2/association/singleTissue{kind}"
 page=0; out=[];receipt=[]
 try:
  while True:
   url=path+'?'+urllib.parse.urlencode({'datasetId':'gtex_v10','variantId':variantid,'itemsPerPage':250,'page':page})
   p=subprocess.run(['curl','-fsSL','--max-time','35',url],capture_output=True,timeout=40)
   if p.returncode:raise ValueError('curl '+str(p.returncode)+' '+p.stderr.decode(errors='replace')[:200])
   data=json.loads(p.stdout);info=data['paging_info'];v=data['data']; key=f"{x['variant']}_{kind}_{ref}_{alt}_p{page}";(raw/(key+'.json')).write_bytes(p.stdout);receipt.append({'url':url,'sha256':hashlib.sha256(p.stdout).hexdigest(),'page':page,'n':len(v),'paging':info})
   for row in v:
    assert row.get('variantId')==variantid and row.get('datasetId')=='gtex_v10'
    if row.get('snpId') not in (None,'',x['variant']):raise ValueError('rsID mismatch '+str(row.get('snpId')))
    out.append({'candidate':x['variant'],'qtl_type':kind,'query_variant':variantid,'row':row})
   count=int(info['numberOfPages']);page+=1
   if page>=count:break
  return {'candidate':x['variant'],'kind':kind,'variantid':variantid,'status':'success','receipts':receipt,'results':out}
 except Exception as e:return {'candidate':x['variant'],'kind':kind,'variantid':variantid,'status':'error','error':str(e),'receipts':receipt,'results':out}
tasks=[(x,k,ref,alt) for x in targets for k in ('Eqtl','Sqtl') for ref,alt in ((x['a'],x['b']),(x['b'],x['a']))]
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
 results=list(ex.map(query,tasks))
pathlib.Path('outputs/stage3_gtex_v10_query_receipts.json').write_text(json.dumps(results,indent=2)+'\n')
rows=[]
for x in results:
 for v in x['results']:
  a=v['row'];rows.append({'candidate':v['candidate'],'qtl_type':v['qtl_type'],'gtex_variant_id':v['query_variant'],'snp_id':a.get('snpId',''),'gene_symbol':a.get('geneSymbol',''),'gencode_id':a.get('gencodeId',''),'tissue':a.get('tissueSiteDetailId',''),'p_value':a.get('pValue',''),'normalized_effect_size':a.get('nes',''),'molecular_phenotype':a.get('phenotypeId',a.get('sQtlId','')),'dataset_id':a.get('datasetId',''),'note':'Precomputed portal association; does not establish causal colocalization or signed alignment to counted allele'})
fields=['candidate','qtl_type','gtex_variant_id','snp_id','gene_symbol','gencode_id','tissue','p_value','normalized_effect_size','molecular_phenotype','dataset_id','note']
with open('outputs/stage3_gtex_v10_candidate_qtls.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
print('queries',len(results),'errors',sum(x['status']!='success' for x in results),'association rows',len(rows))
for t in targets:
 mine=[r for r in results if r['candidate']==t['variant']]; print(t['variant'],'e',sum(len(x['results']) for x in mine if x['kind']=='Eqtl'),'s',sum(len(x['results']) for x in mine if x['kind']=='Sqtl'),'errors',sum(x['status']!='success' for x in mine))

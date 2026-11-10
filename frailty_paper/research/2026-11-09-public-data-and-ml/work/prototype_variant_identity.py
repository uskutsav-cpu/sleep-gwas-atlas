import csv,json,subprocess,concurrent.futures,hashlib
rows=list(csv.DictReader(open('outputs/unvalidated_15_variant_score_prototype.csv')))
def get(x):
 v=x['variant']; url='https://grch37.rest.ensembl.org/variation/human/'+v+'?content-type=application/json'; o={'variant':v,'chr_expected':x['chr'],'pos_expected':x['position_GRCh37_inferred'],'A1':x['counted_allele'],'A2':x['other_allele'],'status':'','alleles':'','coordinate_match':'','alleles_in_reference':'','mapping_count':'','source_url':url,'response_sha256':''}
 try:
  b=subprocess.run(['curl','--fail','--silent','--show-error','--max-time','30',url],capture_output=True,check=True).stdout; d=json.loads(b);maps=[m for m in d.get('mappings',[]) if m.get('assembly_name')=='GRCh37' and str(m.get('seq_region_name'))==x['chr'] and str(m.get('start'))==x['position_GRCh37_inferred']];o['response_sha256']=hashlib.sha256(b).hexdigest();o['mapping_count']=len(maps)
  if len(maps)==1:
   m=maps[0]; aa=m.get('allele_string','').split('/');o.update(status='matched_reference_position',alleles=m.get('allele_string',''),coordinate_match='yes',alleles_in_reference=str(x['counted_allele'] in aa and x['other_allele'] in aa),n_alleles=len(aa))
  else:o.update(status='no_unique_reference_mapping',coordinate_match='no')
 except Exception as e:o['status']='request_failed'
 return o
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:out=list(ex.map(get,rows))
fields=['variant','chr_expected','pos_expected','A1','A2','status','alleles','n_alleles','coordinate_match','alleles_in_reference','mapping_count','source_url','response_sha256']
with open('outputs/prototype_variant_identity_reference.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(out)
print('matched',sum(x['status']=='matched_reference_position' for x in out),'multi-allelic',sum(x.get('n_alleles',0)>2 for x in out),'mismatch',sum(x['alleles_in_reference']=='False' for x in out))
for x in out:
 if x.get('n_alleles',0)>2 or x['status']!='matched_reference_position' or x['alleles_in_reference']=='False':print(x['variant'],x['status'],x['alleles'],x['alleles_in_reference'])

import csv,json,subprocess,concurrent.futures,hashlib,pathlib
markers=list(csv.DictReader(open('outputs/unvalidated_15_variant_score_prototype.csv')))
def get(x):
 v=x['variant'];u='https://rest.ensembl.org/variation/human/'+v+'?content-type=application/json'
 out={'variant':v,'candidate_chr':x['chr'],'candidate_pos_grch37':x['position_GRCh37_inferred'],'Ensembl_GRCh38_mappings':'','status':'','url':u,'sha256':''}
 try:
  b=subprocess.run(['curl','--fail','--silent','--show-error','--max-time','30',u],capture_output=True,check=True).stdout;d=json.loads(b);maps=[{'chr':m['seq_region_name'],'pos':m['start'],'alleles':m.get('allele_string'),'strand':m.get('strand'),'assembly':m.get('assembly_name')} for m in d.get('mappings',[]) if str(m.get('seq_region_name'))==x['chr'] and str(m.get('assembly_name','')).startswith('GRCh38')];out.update(status='success',Ensembl_GRCh38_mappings=json.dumps(maps,separators=(',',':')),sha256=hashlib.sha256(b).hexdigest())
 except Exception:out['status']='request_failed'
 return out
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as e:r=list(e.map(get,markers))
with open('work/variant_reference_grch38.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=r[0]);w.writeheader();w.writerows(r)
print('success',sum(x['status']=='success' for x in r),'with one mapping',sum(len(json.loads(x['Ensembl_GRCh38_mappings']))==1 for x in r if x['status']=='success'))

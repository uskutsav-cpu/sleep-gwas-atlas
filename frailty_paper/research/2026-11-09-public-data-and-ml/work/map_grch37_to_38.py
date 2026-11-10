import csv,json,subprocess,concurrent.futures,hashlib
markers=list(csv.DictReader(open('outputs/unvalidated_15_variant_score_prototype.csv')))
def get(x):
 v=x['variant'];chr=x['chr'];pos=x['position_GRCh37_inferred'];u=f'https://grch37.rest.ensembl.org/map/human/GRCh37/{chr}:{pos}..{pos}/GRCh38?content-type=application/json'
 out={'variant':v,'chr_grch37':chr,'pos_grch37':pos,'grch38_chr':'','grch38_pos':'','mapping_strand':'','status':'','source_url':u,'response_sha256':''}
 try:
  b=subprocess.run(['curl','--fail','--silent','--show-error','--max-time','20',u],capture_output=True,check=True).stdout;d=json.loads(b);m=d.get('mappings',[]);out['response_sha256']=hashlib.sha256(b).hexdigest()
  if len(m)==1:
   a=m[0]['mapped']; out.update(status='unique_mapping',grch38_chr=a['seq_region_name'],grch38_pos=a['start'],mapping_strand=a['strand'])
  else:out['status']='not_unique_or_unmapped'
 except Exception:out['status']='request_failed'
 return out
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as e:r=list(e.map(get,markers))
with open('work/grch37_to_grch38_mapping.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=r[0]);w.writeheader();w.writerows(r)
print('unique mapping',sum(x['status']=='unique_mapping' for x in r),'non-success',[(x['variant'],x['status']) for x in r if x['status']!='unique_mapping'])

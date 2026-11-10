import csv,json,subprocess,hashlib,pathlib,concurrent.futures
m=list(csv.DictReader(open('outputs/shared22_eight_source_evidence_matrix.csv')));ex=next(x for x in csv.DictReader(open('outputs/integrated_15_candidate_review.csv')) if x['variant']=='rs323509');tg=[{'variant':x['variant'],'chr':x['chr_GRCh38'],'pos':x['position_GRCh38']} for x in m]+[{'variant':'rs323509','chr':ex['GRCh37_chr'],'pos':ex['GRCh38_pos']}];raw=pathlib.Path('work/stage3_ensembl_regulatory_raw');raw.mkdir(exist_ok=True)
def run(x):
 u=f"https://rest.ensembl.org/overlap/region/human/{x['chr']}:{x['pos']}-{x['pos']}?feature=regulatory;content-type=application/json";c=subprocess.run(['curl','-fsSL','--max-time','40',u],capture_output=True)
 if c.returncode:return {**x,'status':'error','error':c.stderr.decode(errors='replace')[:100],'url':u}
 try:
  data=json.loads(c.stdout); assert isinstance(data,list)
  for z in data:assert int(z['start'])<=int(x['pos'])<=int(z['end']) and z['seq_region_name']==x['chr']
  (raw/(x['variant']+'.json')).write_bytes(c.stdout)
  return {**x,'status':'success','url':u,'sha256':hashlib.sha256(c.stdout).hexdigest(),'features':data}
 except Exception as e:return {**x,'status':'error','url':u,'error':str(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:r=list(pool.map(run,tg))
pathlib.Path('outputs/stage3_ensembl_regulatory_overlap_receipts.json').write_text(json.dumps({'release':json.load(open('work/stage3_ensembl_release.json')),'queries':r},indent=2)+'\n')
out=[]
for x in r:
 for z in x.get('features',[]):out.append({'variant':x['variant'],'chr_GRCh38':x['chr'],'position_GRCh38':x['pos'],'regulatory_feature_id':z.get('id',''),'feature_type':z.get('feature_type',''),'description':z.get('description',''),'start':z.get('start',''),'end':z.get('end',''),'source_url':x['url'],'evidence_status':'Positional regulatory-feature overlap only; not tissue-specific causal function or a linked effector gene'})
cols=['variant','chr_GRCh38','position_GRCh38','regulatory_feature_id','feature_type','description','start','end','source_url','evidence_status']
with open('outputs/stage3_ensembl_regulatory_overlap.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(out)
print('queries',len(r),'success',sum(x['status']=='success' for x in r),'errors',[(x['variant'],x.get('error')) for x in r if x['status']!='success'],'variants_with_overlap',len(set(x['variant'] for x in out)),'feature rows',len(out))

import csv,json,pathlib,hashlib,math
p=pathlib.Path(__file__).resolve().parent;result=[];files=[]
def vec(obj,key):
 x=obj[key];return [float(v[0] if isinstance(v,list) else v) for v in x]
def se(v):
 n=len(v);m=sum(v)/n;return math.sqrt((n-1)/n*sum((x-m)**2 for x in v))
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
for sleep in ('r13','ukb','mvp'):
 d={}
 for name in ('old_FI','old_FFS','mak_FI','mak_FP'):
  f=p/f'mak_{sleep}_{name}.exact.json';files.append(f);r=json.loads(f.read_text());assert len(r)==1;d[name]=r[0]
 x=list(d.values());assert len({q['row_order_SNP_sha256'] for q in x})==1 and len({q['n_snps'] for q in x})==1
 ref=x[0]['hsq2'];assert all(q['hsq2']['separators']==ref['separators'] and q['hsq1']['separators']==ref['separators'] and len(q['delete_rg'])==200 for q in x)
 h=vec(ref,'delete_total');assert all(max(abs(a-b) for a,b in zip(h,vec(q['hsq2'],'delete_total')))<1e-12 and abs(q['hsq2']['total']-ref['total'])<1e-12 for q in x)
 for name,q in d.items():assert abs(se(vec(q,'delete_rg'))-q['rg_se'])<1e-10 and all(math.isfinite(v) for v in vec(q,'delete_rg'))
 definitions=[('Mak FP − Mak FI',{'mak_FP':1,'mak_FI':-1}),('Ye FFS − Atkins FI',{'old_FFS':1,'old_FI':-1}),('Mak contrast − original contrast',{'mak_FP':1,'mak_FI':-1,'old_FFS':-1,'old_FI':1})]
 for label,coef in definitions:
  delta=sum(c*d[k]['rg'] for k,c in coef.items());v=[sum(c*vec(d[k],'delete_rg')[i] for k,c in coef.items()) for i in range(200)];err=se(v);z=delta/err;pv=math.erfc(abs(z)/math.sqrt(2));family=3 if len(coef)==4 else 6
  result.append({'sleep_source':sleep,'contrast':label,'n_snps':x[0]['n_snps'],'difference':delta,'paired_se':err,'ci95_low':delta-1.95996398454*err,'ci95_high':delta+1.95996398454*err,'p_unadjusted':pv,'bonferroni_p':min(1,pv*family),'correction_family':family})
with (p/'audited_contrasts.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(result[0]),delimiter='\t');w.writeheader();w.writerows(result)
receipt={'source_json_sha256':{str(f):sha(f) for f in files},'n_source_fits':len(files),'families':'six within-pair contrasts and three difference-of-differences; post-hoc and not exhaustive','checks':'same row-order SNP hash, counts, block boundaries, shared insomnia heritability and delete values, 200-block SEs','limitation':'Arithmetic audit of saved native records; original source QC and independence not established.'};(p/'audit_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print('verified',len(files),'fit records,',len(result),'contrasts')

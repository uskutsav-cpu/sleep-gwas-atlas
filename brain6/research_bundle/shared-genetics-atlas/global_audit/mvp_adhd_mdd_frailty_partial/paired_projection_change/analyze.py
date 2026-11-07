import csv,hashlib,json,math,re
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent;atlas=next(q for q in p.parents if q.name=='shared-genetics-atlas');base=atlas/'multiscale_embedding/hfrs_extension/heldout_sleep/group_holdout/bounded_model/mvp_feature_transfer/common_support_run';dirs=[base/'jackknife_uncertainty',atlas/'global_audit/mvp_adhd_frailty_partial',p.parent];logs=[next(d.glob('ldsc-*.log')) for d in dirs];texts=[x.read_text() for x in logs];assert all('Analysis finished' in s for s in texts) and [s.count('991707 SNPs with valid alleles.') for s in texts]==[9,3,4]
pat=r'Genetic Correlation: ([\d.\-]+) \(([\d.]+)\)\nZ-score: [\d.\-]+\nP: ([\deE.\-]+)';matches=[re.findall(pat,s) for s in texts];assert [len(z) for z in matches]==[9,3,4];names=['FI','FFS','HFRS'];source_names=[['adhd','mdd','scz','bipolar','alz','parkinson']+names,names,['adhd']+names];points=[{n:(float(x[0]),float(x[1])) for n,x in zip(ns,ma)} for ns,ma in zip(source_names,matches)];files=logs+[base/'input_verification.json',dirs[1]/'results.tsv',dirs[2]/'results.tsv'];r=json.loads((base/'input_verification.json').read_text());assert len(r['files'])==10 and all(x['rows']==991707 for x in r['files'].values())

def read(path):
 with path.open(newline='') as f:return {r['frailty']:r for r in csv.DictReader(f,delimiter='\t')}
one=read(dirs[1]/'results.tsv');two=read(dirs[2]/'results.tsv');seen_h={}
def vec(d,first,second,pt):
 vals=[]
 for typ in ('hsq1','hsq2','gencov'):
  f=list(d.glob('*'+first+'.sumstats.gz_'+second+'.sumstats.gz.'+typ+'.delete'));assert len(f)==1;files.append(f[0]);x=np.loadtxt(f[0]);assert x.shape==(200,) and np.isfinite(x).all();vals.append(x)
 h1,h2,cov=vals;assert (h1>0).all() and (h2>0).all()
 for n,h in ((first,h1),(second,h2)):
  if n in seen_h:assert np.array_equal(seen_h[n],h),n+' h2 vectors differ'
  else:seen_h[n]=h
 v=cov/np.sqrt(h1*h2);se=float(np.sqrt(199/200*np.sum((v-v.mean())**2)));assert abs(se-pt[1])<0.00006,(first,second,se,pt);return v
ia=vec(dirs[0],'MVP','adhd',points[0]['adhd']);im=vec(dirs[0],'MVP','mdd',points[0]['mdd']);am=vec(dirs[2],'mdd','adhd',points[2]['adhd']);out=[];minimum=1
for n in names:
 inf=vec(dirs[0],'MVP',n,points[0][n]);af=vec(dirs[1],'adhd',n,points[1][n]);mf=vec(dirs[2],'mdd',n,points[2][n]);full=[points[0]['adhd'][0],points[0]['mdd'][0],points[2]['adhd'][0],points[0][n][0],points[1][n][0],points[2][n][0]]
 def calc(vals):
  global minimum
  ia,im,am,inf,af,mf=vals;M=np.array([[1,ia,im,inf],[ia,1,am,af],[im,am,1,mf],[inf,af,mf,1]],float);eig=float(np.linalg.eigvalsh(M).min());minimum=min(minimum,eig);assert eig>1e-8
  cond=M[1:3,1:3];u=M[0,1:3];v=M[3,1:3];ru=1-u@np.linalg.solve(cond,u);rv=1-v@np.linalg.solve(cond,v);assert ru>0 and rv>0;z=(inf-u@np.linalg.solve(cond,v))/np.sqrt(ru*rv);assert abs(z)<=1;return float(z)
 estimate=calc(full);deletes=np.array([calc(x) for x in zip(ia,im,am,inf,af,mf)]);se2=float(np.sqrt(199/200*np.sum((deletes-deletes.mean())**2)))
 def single(x,y,z):return (z-x*y)/np.sqrt((1-x*x)*(1-y*y))
 onepoint=float(single(full[0],full[4],full[3]));onedel=single(ia,af,inf);se1=float(np.sqrt(199/200*np.sum((onedel-onedel.mean())**2)))
 assert abs(onepoint-float(one[n]['partial_genetic_rg']))<1e-10 and abs(se1-float(one[n]['paired_jackknife_se']))<1e-10
 assert abs(estimate-float(two[n]['partial_rg_given_adhd_mdd']))<1e-10 and abs(se2-float(two[n]['paired_jackknife_se']))<1e-10
 d=estimate-onepoint;dd=deletes-onedel;se=float(np.sqrt(199/200*np.sum((dd-dd.mean())**2)));pv=math.erfc(abs(d/se)/math.sqrt(2));out.append({'frailty':n,'adhd_only':onepoint,'adhd_and_mdd':estimate,'difference_two_minus_one':d,'paired_se':se,'normal_95_low':d-1.96*se,'normal_95_high':d+1.96*se,'unadjusted_p':pv,'bonferroni_p_three':min(1,3*pv)})
with (p/'results.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=out[0],delimiter='\t');w.writeheader();w.writerows(out)
(p/'receipt.json').write_text(json.dumps({'snps':991707,'blocks':200,'min_eigenvalue':minimum,'h2_delete_vectors_match':sorted(seen_h),'input_sha256':{str(x.relative_to(atlas)):hashlib.sha256(x.read_bytes()).hexdigest() for x in files},'scope':'post hoc genetic projections not individual diagnoses, causality or interaction'},indent=2)+'\n');print(json.dumps(out,indent=2));print('min_eigenvalue',minimum)

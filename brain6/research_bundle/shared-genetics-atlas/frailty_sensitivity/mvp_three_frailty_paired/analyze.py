import csv,hashlib,json,math,re
from pathlib import Path
import numpy as np
from scipy.special import ndtr
p=Path(__file__).resolve().parent;atlas=next(q for q in p.parents if q.name=='shared-genetics-atlas');base=atlas/'multiscale_embedding/hfrs_extension/heldout_sleep/group_holdout/bounded_model/mvp_feature_transfer/common_support_run';j=base/'jackknife_uncertainty';logs=list(j.glob('ldsc-*.log'));assert len(logs)==1;log=logs[0];s=log.read_text();assert 'Analysis finished' in s and s.count('991707 SNPs with valid alleles.')==9
matches=re.findall(r'Genetic Correlation: ([\d.\-]+) \(([\d.]+)\)\nZ-score: [\d.\-]+\nP: ([\deE.\-]+)',s);assert len(matches)==9
names=['adhd','mdd','scz','bipolar','alz','parkinson','FI','FFS','HFRS'];points={n:float(m[0]) for n,m in zip(names,matches)};ses={n:float(m[1]) for n,m in zip(names,matches)};arr={};first=[];inputs=[log,base/'input_verification.json']
for n in ('FI','FFS','HFRS'):
 parts=[]
 for part in ('hsq1','hsq2','gencov'):
  files=list(j.glob('*MVP.sumstats.gz_'+n+'.sumstats.gz.'+part+'.delete'));assert len(files)==1;inputs.append(files[0]);v=np.loadtxt(files[0]);assert v.shape==(200,) and np.isfinite(v).all();parts.append(v)
 h1,h2,cov=parts;assert (h1>0).all() and (h2>0).all();first.append(h1);v=cov/np.sqrt(h1*h2);se=float(np.sqrt(199/200*np.sum((v-v.mean())**2)));assert abs(se-ses[n])<0.00006;arr[n]=v
assert all(np.array_equal(first[0],x) for x in first[1:]);verify=json.loads((base/'input_verification.json').read_text());assert len(verify['files'])==10 and all(x['rows']==991707 for x in verify['files'].values());out=[]
for left,right in (('FI','FFS'),('FI','HFRS'),('FFS','HFRS')):
 d=points[left]-points[right];v=arr[left]-arr[right];se=float(np.sqrt(199/200*np.sum((v-v.mean())**2)));assert se>0;z=d/se;pv=float(2*ndtr(-abs(z)));out.append({'contrast':left+' minus '+right,'difference':d,'paired_jackknife_se':se,'normal_95_low':d-1.96*se,'normal_95_high':d+1.96*se,'unadjusted_p':pv,'bonferroni_p_three':min(1,3*pv)})
with (p/'contrasts.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=out[0],delimiter='\t');w.writeheader();w.writerows(out)
(p/'receipt.json').write_text(json.dumps({'matched_snps':991707,'blocks':200,'inputs_sha256':{str(x.relative_to(atlas)):hashlib.sha256(x.read_bytes()).hexdigest() for x in inputs},'identical_mvp_h2_delete_vectors':True,'point_estimates':'rounded LDSC log values','scope':'post-hoc three-contrast family, not project-wide'},indent=2)+'\n');print(json.dumps(out,indent=2))

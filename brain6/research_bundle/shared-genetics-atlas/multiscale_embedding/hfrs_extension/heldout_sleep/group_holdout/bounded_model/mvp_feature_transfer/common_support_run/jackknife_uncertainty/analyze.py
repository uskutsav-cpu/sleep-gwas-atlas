import csv,glob,hashlib,json,re
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent; parent=p.parent;atlas=next(q for q in p.parents if q.name=='shared-genetics-atlas');source=atlas/'multiscale_embedding/hfrs_extension/source_matrix.tsv';pred=parent/'predictions.tsv';features=parent/'features.tsv';verified=parent/'input_verification.json';assert json.loads(verified.read_text())['result'].startswith('all files read')
logs=list(p.glob('ldsc-*.log'));assert len(logs)==1; log=logs[0];s=log.read_text();assert 'Analysis finished' in s and s.count('991707 SNPs with valid alleles.')==9
names=['adhd','mdd','scz','bipolar','alz','parkinson','FI','FFS','HFRS']; matches=re.findall(r'Genetic Correlation: ([\d.\-]+) \(([\d.]+)\)\nZ-score: [\d.\-]+\nP: ([\deE.\-]+)',s);assert len(matches)==9
vectors=[];hsqfirst=[];files=[];ses=[]
for i,name in enumerate(names):
 parts=[]
 for part in ('hsq1','hsq2','gencov'):
  found=list(p.glob('*MVP.sumstats.gz_'+name+'.sumstats.gz.'+part+'.delete'));assert len(found)==1; files.append(found[0]);v=np.loadtxt(found[0]);assert v.shape==(200,) and np.isfinite(v).all();parts.append(v)
 h1,h2,cov=parts;assert (h1>0).all() and (h2>0).all();hsqfirst.append(h1);v=cov/np.sqrt(h1*h2);assert np.isfinite(v).all();err=float(np.sqrt(199/200*np.sum((v-v.mean())**2)));assert abs(err-float(matches[i][1]))<0.00006,(name,err,matches[i][1]);vectors.append(v);ses.append(err)
assert all(np.array_equal(hsqfirst[0],x) for x in hsqfirst[1:]); V=np.stack(vectors,axis=1)
def read(path):
 with path.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
r=read(source);olds=read(pred);fs=read(features); assert [x['feature'] for x in fs]==names[:6] and [x['target'] for x in olds]==names[6:]
assert all(abs(float(matches[i][0])-float(fs[i]['rg']))<1e-9 for i in range(6))
assert all(abs(float(matches[6+i][0])-float(olds[i]['observed_mvp_rg']))<1e-9 for i in range(3))
lookup={(x['sleep_trait'],x['outcome']):float(x['rg']) for x in r};train=[n for n in dict.fromkeys(x['sleep_trait'] for x in r) if n!='insomnia'];X=np.array([[lookup[n,k] for k in names[:6]] for n in train]);Y=np.array([[lookup[n,k] for k in names[6:]] for n in train]);mu=X.mean(0);sd=X.std(0);sd=np.where(sd>1e-12,sd,1);A=(X-mu)/sd;point=np.array([float(x['rg']) for x in fs]);out=[]
for j,t in enumerate(names[6:]):
 lam=float(olds[j]['lambda']);z=np.arctanh(Y[:,j]);b=np.linalg.lstsq(np.vstack((A,np.sqrt(lam)*np.eye(6))),np.r_[z-z.mean(),np.zeros(6)],rcond=None)[0];fn=lambda x:np.tanh(z.mean()+((x-mu)/sd)@b);estimate=float(fn(point));assert abs(estimate-float(olds[j]['mvp_feature_bounded_prediction']))<1e-10
 residual=estimate-float(olds[j]['observed_mvp_rg']); deletes=fn(V[:,:6])-V[:,6+j];se=float(np.sqrt(199/200*np.sum((deletes-deletes.mean())**2)));out.append({'target':t,'point_residual_prediction_minus_observed':residual,'conditional_jackknife_se':se,'conditional_normal_low':residual-1.96*se,'conditional_normal_high':residual+1.96*se})
with (p/'residual_uncertainty.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=out[0],delimiter='\t');w.writeheader();w.writerows(out)
(p/'receipt.json').write_text(json.dumps({'delete_blocks':200,'feature_and_target_rg_jackknife_se':dict(zip(names,ses)),'identical_first_trait_h2_delete_vectors':True,'input_hashes':{str(x.relative_to(atlas)):hashlib.sha256(x.read_bytes()).hexdigest() for x in [source,pred,features,verified,log,*files]},'qualification':'conditional on fixed model and LDSC data; not total predictive uncertainty'},indent=2)+'\n');print(json.dumps(out,indent=2))

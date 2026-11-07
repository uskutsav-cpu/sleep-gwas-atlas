import csv,hashlib,json,re
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent
atlas=next(q for q in p.parents if q.name=='shared-genetics-atlas')
source=atlas/'multiscale_embedding/hfrs_extension/source_matrix.tsv'; saved=p.parent.parent/'predictions.tsv'; targetlog=None
logs=list(p.glob('ldsc-*.log'));receipts=list(p.glob('mvp_stream_receipt.json')); assert len(logs)==len(receipts)==1
log=logs[0];targetlog=log;receipt=receipts[0];prior=atlas/'replication/mvp_finngen_genomewide/stream_receipt.json';assert json.loads(receipt.read_text())==json.loads(prior.read_text())
s=log.read_text();assert 'Analysis finished' in s
features=['adhd','mdd','scz','bipolar','alz','parkinson'];targets=['FI','FFS','HFRS']
positions=[s.index('/'+x+'.sumstats.gz ...') for x in features];assert positions==sorted(positions)
allmatches=re.findall(r'Genetic Correlation: ([\d.\-]+) \(([\d.]+)\)\nZ-score: [\d.\-]+\nP: ([\deE.\-]+)',s);assert len(allmatches)==9;matches=allmatches[:6]
v=np.array([float(m[0]) for m in matches]);assert np.isfinite(v).all() and np.all(abs(v)<1)
assert s.count('991707 SNPs with valid alleles.')==9
outcomes=allmatches[6:]
with source.open(newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
lookup={(r['sleep_trait'],r['outcome']):float(r['rg']) for r in rows}; names=list(dict.fromkeys(r['sleep_trait'] for r in rows)); train=[n for n in names if n!='insomnia'];assert len(train)==11
X=np.array([[lookup[n,k] for k in features] for n in train]);Y=np.array([[lookup[n,k] for k in targets] for n in train]);old=np.array([lookup['insomnia',k] for k in features]);mu=X.mean(0);sd=X.std(0);sd=np.where(sd>1e-12,sd,1);A=(X-mu)/sd
with saved.open(newline='') as f:prev={r['target']:r for r in csv.DictReader(f,delimiter='\t') if r['sleep_trait']=='insomnia'}
assert set(prev)==set(targets)
def knn(x):return np.argsort(((A-(x-mu)/sd)**2).sum(1),kind='stable')[:3]
out=[]
for j,t in enumerate(targets):
 r=prev[t];lam=float(r['lambda']); z=np.arctanh(Y[:,j]); zm=z.mean();coef=np.linalg.lstsq(np.vstack([A,np.sqrt(lam)*np.eye(6)]),np.r_[z-zm,np.zeros(6)],rcond=None)[0]
 predict=lambda x: float(np.tanh(zm+((x-mu)/sd)@coef))
 assert abs(predict(old)-float(r['bounded_ridge']))<1e-10
 assert abs(float(Y[:,j].mean())-float(r['training_mean']))<1e-10
 assert abs(float(Y[knn(old),j].mean())-float(r['knn3']))<1e-10
 observed=float(outcomes[j][0]); prediction=predict(v); mean=float(Y[:,j].mean()); nearest=knn(v); near=float(Y[nearest,j].mean())
 out.append({'target':t,'observed_mvp_rg':observed,'observed_mvp_se':float(outcomes[j][1]),'mvp_feature_bounded_prediction':prediction,'bounded_abs_error':abs(prediction-observed),'mean_prediction':mean,'mean_abs_error':abs(mean-observed),'knn3_prediction':near,'knn3_abs_error':abs(near-observed),'knn3_neighbors':','.join(train[i] for i in nearest),'questionnaire_proxy_prediction':float(r['bounded_ridge']),'lambda':lam})
with (p/'predictions.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=out[0],delimiter='\t');w.writeheader();w.writerows(out)
feature_rows=[{'feature':k,'rg':m[0],'se':m[1],'p':m[2]} for k,m in zip(features,matches)]
with (p/'features.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=feature_rows[0],delimiter='\t');w.writeheader();w.writerows(feature_rows)
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
(p/'prediction_receipt.json').write_text(json.dumps({'input_sha256':{str(f.relative_to(atlas)):sha(f) for f in (source,saved,log,receipt,prior,p/'prepared_receipt.json')},'training_rows':train,'excluded':'insomnia','validation':'Old questionnaire predictions and baselines reproduced to 1e-10 using independent augmented least squares','feature_order':features,'test_interpretation':'post-hoc common-marker sensitivity; not independent validation'},indent=2)+'\n')
print(json.dumps(out,indent=2))

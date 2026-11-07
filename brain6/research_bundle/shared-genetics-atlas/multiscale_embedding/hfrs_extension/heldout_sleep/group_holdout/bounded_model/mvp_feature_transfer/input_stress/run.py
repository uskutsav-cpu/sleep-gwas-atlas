import csv,hashlib,itertools,json
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent; base=p.parent;atlas=next(q for q in p.parents if q.name=='shared-genetics-atlas')
source=atlas/'multiscale_embedding/hfrs_extension/source_matrix.tsv'; features=base/'features.tsv'; predictions=base/'predictions.tsv'
def read(x):
 with x.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
r=read(source); feat=read(features); prev=read(predictions); names=list(dict.fromkeys(x['sleep_trait'] for x in r));train=[x for x in names if x!='insomnia']; keys=['adhd','mdd','scz','bipolar','alz','parkinson'];targets=['FI','FFS','HFRS'];assert [x['feature'] for x in feat]==keys and len(train)==11 and [x['target'] for x in prev]==targets
lookup={(x['sleep_trait'],x['outcome']):float(x['rg']) for x in r};X=np.array([[lookup[n,k] for k in keys] for n in train]);Y=np.array([[lookup[n,k] for k in targets] for n in train]);v=np.array([float(x['rg']) for x in feat]);se=np.array([float(x['se']) for x in feat]);assert np.isfinite(se).all() and np.all(se>0)
mu=X.mean(0);sd=X.std(0);sd=np.where(sd>1e-12,sd,1);a=(X-mu)/sd;out=[];maxdiff=0
for j,t in enumerate(targets):
 prior=prev[j];lam=float(prior['lambda']); z=np.arctanh(Y[:,j]);beta=np.linalg.solve(a.T@a+lam*np.eye(6),a.T@(z-z.mean()));c=float(z.mean()+((v-mu)/sd)@beta);rad=float(np.sum(abs(beta*se/sd))); assert rad>0 and abs(np.tanh(c)-float(prior['mvp_feature_bounded_prediction']))<1e-10
 lo,hi=np.tanh([c-rad,c+rad]); corners=[float(np.tanh(z.mean()+((v+np.array(s)*se-mu)/sd)@beta)) for s in itertools.product((-1,1),repeat=6)];maxdiff=max(maxdiff,abs(min(corners)-lo),abs(max(corners)-hi)); target=float(prior['observed_mvp_rg']);delta=float(np.arctanh(target)-c);mult=abs(delta)/rad; chosen=v+np.sign(delta)*np.sign(beta)*se*mult;check=float(np.tanh(z.mean()+((chosen-mu)/sd)@beta));assert abs(check-target)<1e-10
 out.append({'target':t,'point_prediction':float(np.tanh(c)),'min_1se_box':float(lo),'max_1se_box':float(hi),'observed_target':target,'observed_in_box':bool(lo<=target<=hi),'min_common_se_multiple_to_target':mult,'required_input_outside_valid_range':bool(np.any(abs(chosen)>1)),'maximum_absolute_corner_difference':maxdiff})
assert maxdiff<1e-12
with (p/'ranges.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=out[0],delimiter='\t');w.writeheader();w.writerows(out)
(p/'receipt.json').write_text(json.dumps({'inputs_sha256':{str(x.relative_to(atlas)):hashlib.sha256(x.read_bytes()).hexdigest() for x in (source,features,predictions)},'corners_per_target':64,'maximum_analytic_corner_difference':maxdiff,'qualification':'fixed-model simultaneous sensitivity boxes, not confidence intervals'},indent=2)+'\n')
print(json.dumps(out,indent=2))

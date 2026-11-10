"""Independent CPU/GPU, zero-contrast and finite-difference checks."""
from pathlib import Path
import json,numpy as np
from contrast_network import mx,initialize,forward,objective
WORK=Path(__file__).resolve().parent;OUT=WORK.parents[1]/'outputs/sleep_frailty_ai_v3'
rng=np.random.default_rng(921);dc,dv=9,7;c=rng.normal(size=(4,dc)).astype('float32');v=np.abs(rng.normal(size=(4,dv))).astype('float32');sid=np.array([0,1,0,1],dtype='int32');pw=np.array([[3,2],[1,1.5]],dtype='float32');p=initialize(dc,dv,42,np.array([.1,.2]),np.array([[.02,.1],[0,.3]]),pw);pm={k:mx.array(a) for k,a in p.items()};log=np.log(pw)
a=forward(pm,mx.array(c),mx.array(v),mx.array(sid),mx.array(log));b=forward(p,c,v,sid,log,backend='numpy');errors={}
for k in a:
    error=float(np.max(np.abs(np.array(a[k])-b[k])));errors[k]=error;assert error<2e-5,(k,error)
null=forward(pm,mx.array(c),mx.zeros(v.shape),mx.array(sid),mx.array(log));assert np.max(np.abs(np.array(null['variant_evidence_logits'])))==0 and np.max(np.abs(np.array(null['absolute_effect'])))==0
assert np.all(np.array(a['conditional_scores'])[sid==1,0]==0)
constants={'log_pw':mx.array(log),'pw':mx.array(pw),'study_weights':mx.array([1.,1.]),'group_weights':mx.array([[1.,1.],[0.,1.]])};masks=[mx.ones((4,128)),mx.ones((4,32)),mx.ones((4,128)),mx.ones((4,32))]
args=[mx.array(c),mx.array(v),mx.array(sid),mx.array([0.,1.,1.,0.]),mx.array([0,1,1,0]),mx.array([0.,.2,.5,-.1]),mx.array([.2,.3,.4,.1]),mx.array([0.,1.,-1.,.5]),mx.ones(4),masks,constants,.05]
loss=lambda par:objective(par,*args);gradient=mx.grad(loss)(pm);checks=[]
for name,index in [('conditional.b',(1,)),('variant.w2',(3,5)),('context.w2',(2,4))]:
    original=p[name][index];eps=.02 if name=='context.w2' else .001;plus={k:mx.array(x) for k,x in p.items()};minus={k:mx.array(x) for k,x in p.items()};x=p[name].copy();x[index]=original+eps;plus[name]=mx.array(x);x=p[name].copy();x[index]=original-eps;minus[name]=mx.array(x);numeric=(float(loss(plus))-float(loss(minus)))/(2*eps);actual=float(np.array(gradient[name])[index]);assert abs(actual-numeric)<max(2e-5,.03*max(abs(actual),abs(numeric))),(name,actual,numeric);checks.append({'parameter':name,'index':index,'autodiff':actual,'finite_difference':numeric,'absolute_difference':abs(actual-numeric)})
receipt={'cpu_gpu_max_errors':errors,'zero_contrast_creates_no_variant_evidence_or_effect_magnitude':True,'measured_activity_not_an_inference_argument':True,'Siraj_inactive_branch_zero_by_label_definition':True,'finite_difference_checks':checks,'scope':'Implementation/structural checks, not biological or clinical validation'}
(OUT/'contrast_implementation_check.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

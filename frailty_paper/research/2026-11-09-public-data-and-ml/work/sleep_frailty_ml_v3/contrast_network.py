"""A context gate and a bias-free variant-contrast network.

Measured activity is never a forward input. At zero DNA contrast, context
cannot create a variant evidence term. A nonzero statistical-call baseline is
still possible: this is not a deterministic mechanistic model or clinical risk.
"""
from pathlib import Path
import sys,numpy as np
OLD=Path(__file__).resolve().parents[1]/'sleep_frailty_ml'
sys.path.insert(0,str(OLD/'runtime'));import mlx.core as mx

def gelu(x,lib):return .5*x*(1+lib.tanh(.7978845608028654*(x+.044715*x**3)))
def sigmoid(x,lib):return 1/(1+lib.exp(-lib.clip(x,-30,30)))

def initialize(dc,dv,seed,gate_prevalence,conditional_prevalence,positive_weight):
    def normal(offset,shape):return np.random.default_rng(seed+offset).normal(0,1/np.sqrt(shape[0]),shape).astype('float32')
    p={'context.w1':normal(1,(dc,128)),'context.b1':np.zeros(128,dtype='float32'),'context.w2':normal(2,(128,32)),'context.b2':np.zeros(32,dtype='float32'),'variant.w1':normal(3,(dv,128)),'variant.w2':normal(4,(128,32)),'modulation.w':normal(5,(32,32))*.1,'modulation.b':np.zeros(32,dtype='float32')}
    p['gate.w']=normal(6,(32,2))*.01;p['gate.b']=np.log(gate_prevalence/(1-gate_prevalence)).astype('float32')
    q=conditional_prevalence;weighted=q*positive_weight/(1-q+q*positive_weight);weighted=np.clip(weighted,1e-6,1-1e-6)
    p['conditional.w']=normal(7,(32,4))*.01;p['conditional.b']=np.log(weighted/(1-weighted)).reshape(-1).astype('float32')
    p['activity.w']=normal(8,(32,2))*.01;p['activity.b']=np.zeros(2,dtype='float32');p['effect.w']=normal(9,(32,2))*.01;p['avi.w']=np.zeros(32,dtype='float32');p['avi.b']=np.array(0.,dtype='float32')
    return p

def forward(p,context,variant,study,log_positive_weight,masks=None,backend='mlx'):
    lib=mx if backend=='mlx' else np
    c=gelu(context@p['context.w1']+p['context.b1'],lib)
    if masks is not None:c=c*masks[0]
    c=gelu(c@p['context.w2']+p['context.b2'],lib)
    if masks is not None:c=c*masks[1]
    v=gelu(variant@p['variant.w1'],lib)
    if masks is not None:v=v*masks[2]
    v=gelu(v@p['variant.w2'],lib)
    if masks is not None:v=v*masks[3]
    evidence=v*(2*sigmoid(c@p['modulation.w']+p['modulation.b'],lib));row=lib.arange(context.shape[0])
    gate_logits=(c@p['gate.w']+p['gate.b'])[row,study];gate=sigmoid(gate_logits,lib)
    variant_terms=(evidence@p['conditional.w']).reshape((-1,2,2))[row,study,:]
    raw_logits=(evidence@p['conditional.w']+p['conditional.b']).reshape((-1,2,2))[row,study,:]
    conditional=sigmoid(raw_logits-log_positive_weight[study],lib)
    # Siraj's predefined label requires source activity. Its inactive call rate
    # is zero by definition. Salomon has genuine calls outside that stratum.
    conditional=lib.stack([lib.where(study==1,0.,conditional[:,0]),conditional[:,1]],axis=1)
    overall=(1-gate)*conditional[:,0]+gate*conditional[:,1]
    return {'overall':overall,'gate_logits':gate_logits,'predicted_active':gate,'raw_conditional_logits':raw_logits,'conditional_scores':conditional,'variant_evidence_logits':variant_terms,'activity':(c@p['activity.w']+p['activity.b'])[row,study],'absolute_effect':lib.abs((evidence@p['effect.w'])[row,study]),'avi':evidence@p['avi.w']+p['avi.b'],'variant_latent':evidence}

def huber(x):return mx.where(mx.abs(x)<=1,.5*x*x,mx.abs(x)-.5)

def objective(p,c,v,sid,y,active,activity,effect,teacher,teacher_mask,masks,constants,avi_weight):
    z=forward(p,c,v,sid,constants['log_pw'],masks);sw=constants['study_weights'][sid];prob=mx.clip(z['overall'],1e-6,1-1e-6)
    call=mx.mean(sw*(-y*mx.log(prob)-(1-y)*mx.log1p(-prob)))
    row=mx.arange(c.shape[0]);raw=z['raw_conditional_logits'][row,active];pw=constants['pw'][sid,active];gw=constants['group_weights'][sid,active]
    conditional=mx.mean(gw*((1-y)*mx.logaddexp(0,raw)+pw*y*mx.logaddexp(0,-raw)))
    gate=mx.mean(sw*((1-active)*mx.logaddexp(0,z['gate_logits'])+active*mx.logaddexp(0,-z['gate_logits'])))
    act=mx.mean(sw*huber(z['activity']-activity));eff=mx.mean(sw*huber(z['absolute_effect']-effect))
    av=mx.sum(sw*teacher_mask*huber(z['avi']-teacher))/mx.maximum(mx.sum(sw*teacher_mask),1.)
    return call+conditional+.25*gate+.1*act+.1*eff+avi_weight*av

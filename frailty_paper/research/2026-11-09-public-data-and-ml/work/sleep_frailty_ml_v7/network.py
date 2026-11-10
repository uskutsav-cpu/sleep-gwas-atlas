"""A quantitative effect network; DNA contrast0 gives predicted effect0."""
import sys
from common import *
sys.path.insert(0,str(OLD/'runtime'))
import mlx.core as mx
def gelu(x,lib):return .5*x*(1+lib.tanh(.7978845608028654*(x+.044715*x**3)))
def initialize(context_dim,variant_dim,seed,avi=False):
    rng=np.random.default_rng(seed)
    def normal(shape,scale=1):return (rng.normal(0,scale/np.sqrt(shape[0]),shape)).astype('float32')
    base_dim=context_dim-3 if avi else context_dim;cw=normal((base_dim,128))
    if avi:cw=np.r_[cw,np.zeros((3,128),dtype='float32')]
    return {'cw1':cw,'cb1':np.zeros(128,dtype='float32'),'cw2':normal((128,64)),'cb2':np.zeros(64,dtype='float32'),'vw1':normal((variant_dim,128)),'vw2':normal((128,64)),'modw':normal((64,64),.1),'modb':np.zeros(64,dtype='float32'),'effect_w':normal((64,2),.02)}
def forward(p,c,v,masks=None,backend='mlx'):
    lib=mx if backend=='mlx' else np;c=gelu(c@p['cw1']+p['cb1'],lib);v=gelu(v@p['vw1'],lib)
    if masks is not None:c=c*masks[0];v=v*masks[2]
    c=gelu(c@p['cw2']+p['cb2'],lib);v=gelu(v@p['vw2'],lib)
    if masks is not None:c=c*masks[1];v=v*masks[3]
    hidden=v*(1+lib.tanh(c@p['modw']+p['modb']));return lib.abs(hidden@p['effect_w'])
def huber(x,lib=mx):
    a=lib.abs(x);return lib.where(a<=1,.5*x*x,a-.5)
def loss(p,c,v,unit,target,weight,masks,nmain,nlocal,ncross,pair_weight,raw_mae=False,scales=None):
    z=forward(p,c,v,masks);row=mx.arange(len(unit));prediction=z[row,unit]
    if raw_mae:
        scale=scales[unit];prediction=mx.expm1(mx.minimum(prediction*scale,10.))/scale;target=mx.expm1(target*scale)/scale
    error=lambda x:mx.abs(x) if raw_mae else huber(x)
    result=mx.mean(weight[:nmain]*error(prediction[:nmain]-target[:nmain]))
    if nlocal:
        q=prediction[nmain:nmain+2*nlocal].reshape(nlocal,2);t=target[nmain:nmain+2*nlocal].reshape(nlocal,2);u=unit[nmain:nmain+2*nlocal].reshape(nlocal,2)[:,0];w=mx.where(u==0,1.,.15);result=result+pair_weight*mx.mean(w*error((q[:,0]-q[:,1])-(t[:,0]-t[:,1])))
    if ncross:
        start=nmain+2*nlocal;q=prediction[start:start+2*ncross].reshape(ncross,2);t=target[start:start+2*ncross].reshape(ncross,2);result=result+pair_weight*mx.mean(error((q[:,0]-q[:,1])-(t[:,0]-t[:,1])))
    return result

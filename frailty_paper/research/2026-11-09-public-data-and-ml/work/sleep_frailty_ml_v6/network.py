"""Small measured-outcome networks with NumPy/MLX-identical forward rules."""
import sys
from common import *
sys.path.insert(0,str(OLD/'runtime'));import mlx.core as mx
def gelu(x,lib):return .5*x*(1+lib.tanh(.7978845608028654*(x+.044715*x**3)))
def initialize(dim,architecture,context_cols,variant_cols,nstudy,prevalence,seed):
    rng=np.random.default_rng(seed)
    def normal(shape,scale=1):return (rng.normal(0,scale/np.sqrt(shape[0]),shape)).astype('float32')
    if architecture=='mlp':
        p={'w1':normal((dim,128)),'b1':np.zeros(128,dtype='float32'),'w2':normal((128,32)),'b2':np.zeros(32,dtype='float32')}
    else:
        p={'cw1':normal((len(context_cols),128)),'cb1':np.zeros(128,dtype='float32'),'cw2':normal((128,32)),'cb2':np.zeros(32,dtype='float32'),'vw1':normal((len(variant_cols),128)),'vw2':normal((128,32)),'modw':normal((32,32),.1),'modb':np.zeros(32,dtype='float32'),'context_class_w':normal((32,),.01)}
    p.update({'class_w':normal((32,),.01),'class_b':np.array(np.log(prevalence/(1-prevalence)),dtype='float32'),'effect_w':normal((32,nstudy),.01),'effect_b':np.zeros(nstudy,dtype='float32'),'activity_w':normal((32,nstudy),.01),'activity_b':np.zeros(nstudy,dtype='float32')})
    return p
def forward(p,x,study,architecture,context_cols,variant_cols,masks=None,backend='mlx'):
    lib=mx if backend=='mlx' else np
    if architecture=='mlp':
        h=gelu(x@p['w1']+p['b1'],lib)
        if masks is not None:h=h*masks[0]
        h=gelu(h@p['w2']+p['b2'],lib)
        if masks is not None:h=h*masks[1]
        c=h;extra=0
    else:
        c=gelu(x[:,context_cols]@p['cw1']+p['cb1'],lib)
        if masks is not None:c=c*masks[0]
        c=gelu(c@p['cw2']+p['cb2'],lib)
        if masks is not None:c=c*masks[1]
        v=gelu(x[:,variant_cols]@p['vw1'],lib)
        if masks is not None:v=v*masks[2]
        v=gelu(v@p['vw2'],lib)
        if masks is not None:v=v*masks[3]
        h=v*(1+lib.tanh(c@p['modw']+p['modb']));extra=c@p['context_class_w']
    row=lib.arange(len(study))
    return h@p['class_w']+p['class_b']+extra,(h@p['effect_w']+p['effect_b'])[row,study],(c@p['activity_w']+p['activity_b'])[row,study]
def branches(cols,slices,ncells):
    M=(slices['local_motif'][1]-slices['local_motif'][0])//3;context_raw=set(range(0,128))|set(range(512,640))|set(range(1024,1108))|set(range(1192,1192+M))|set(range(1192+2*M,1192+3*M));c=np.array([i for i,k in enumerate(cols) if k in context_raw]+list(range(len(cols),len(cols)+ncells)),dtype='int32');v=np.array([i for i,k in enumerate(cols) if k not in context_raw],dtype='int32');assert len(c)+len(v)==len(cols)+ncells
    return c,v

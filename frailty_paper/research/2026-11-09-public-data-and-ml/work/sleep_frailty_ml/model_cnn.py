from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent/'runtime'))
import mlx.core as mx

def initialize():
    # MLX convolutions use output-channel, kernel, input-channel weight order.
    dims=[(32,15,4),(32,5,32),(64,5,32)]
    p={}
    for i,shape in enumerate(dims):
        p[f'conv{i}.w']=mx.random.normal(shape)*((2/(shape[1]*shape[2]))**.5);p[f'conv{i}.b']=mx.zeros(shape[0])
    p.update({'head.w1':mx.random.normal((192,32))*.04,'head.b1':mx.zeros(32),'head.w2':mx.random.normal((32,))*.04,'head.b2':mx.array(-3.5)})
    return p

def encode(t,p):
    h=mx.eye(4)[t-7]
    for i,k in enumerate([15,5,5]):h=mx.maximum(mx.conv1d(h,p[f'conv{i}.w'],padding=k//2)+p[f'conv{i}.b'],0)
    return h

def pair(ref,alt,p):
    B=ref.shape[0];h=encode(mx.concatenate([ref,alt],axis=0),p);a,b=h[:B],h[B:];diff=mx.abs(a-b)
    return mx.concatenate([mx.mean((a+b)/2,axis=1),mx.mean(diff,axis=1),mx.max(diff,axis=1)],axis=1)

def logits(ref,alt,p):
    f=.5*(pair(ref,alt,p)+pair(17-ref[:,::-1],17-alt[:,::-1],p))
    # Per-example RMS keeps the head away from numerical saturation without
    # sharing statistics between training and held-out examples.
    f=f*mx.rsqrt(mx.mean(f*f,axis=-1,keepdims=True)+1e-6)
    return mx.tanh(f@p['head.w1']+p['head.b1'])@p['head.w2']+p['head.b2']

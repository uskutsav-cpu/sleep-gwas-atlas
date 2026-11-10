"""Compact HyenaDNA port using MLX, checked against model.py before training."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent/'runtime'))
import mlx.core as mx

def linear(x,w,b=None):
    y=x@w.T
    return y if b is None else y+b

def norm(x,w,b):
    m=mx.mean(x,axis=-1,keepdims=True)
    return (x-m)*mx.rsqrt(mx.mean((x-m)**2,axis=-1,keepdims=True)+1e-5)*w+b

def gelu(x):return .5*x*(1+mx.tanh(.7978845608028654*(x+.044715*x**3)))

def kernel(w,p,L):
    z=w[p+'filter_fn.pos_emb.z'][:,:L,:];t=w[p+'filter_fn.pos_emb.t'][:,:L,:]
    f=p+'filter_fn.implicit_filter.';freq=w[f+'1.freq']
    h=mx.sin(linear(z,w[f+'0.weight'],w[f+'0.bias'])*freq)
    for i in [2,4]:h=mx.sin(linear(h,w[f+f'{i}.weight'],w[f+f'{i}.bias'])*freq)
    return linear(h,w[f+'6.weight'])*(mx.exp(-t*mx.abs(w[p+'filter_fn.modulation.deltas']))+.05)

def fftconv(u,k,d):
    L=u.shape[-1]
    # Full complex FFT avoids a failed real-FFT autodiff identity in MLX 0.22.1.
    return mx.real(mx.fft.ifft(mx.fft.fft(u,2*L)*mx.fft.fft(k,2*L)))[...,:L]+u*d[:,None]

def block(x,w,i):
    p=f'layers.{i}.';m=p+'mixer.'
    h=norm(x,w[p+'norm1.weight'],w[p+'norm1.bias'])
    h=linear(h,w[m+'in_proj.weight'],w[m+'in_proj.bias'])
    B,L,D=h.shape;pad=mx.concatenate([mx.zeros((B,2,D),dtype=h.dtype),h],axis=1)
    filt=w[m+'short_filter.weight'][:,0,:]
    uc=sum(pad[:,j:j+L,:]*filt[:,j] for j in range(3))+w[m+'short_filter.bias']
    x0,x1,v=mx.split(uc,3,axis=-1);k=kernel(w,m,L)[0].T
    y=mx.transpose(fftconv(mx.transpose(v*x1,(0,2,1)),k,w[m+'filter_fn.bias']),(0,2,1))*x0
    x=x+linear(y,w[m+'out_proj.weight'],w[m+'out_proj.bias'])
    h=norm(x,w[p+'norm2.weight'],w[p+'norm2.bias'])
    return x+linear(gelu(linear(h,w[p+'mlp.fc1.weight'],w[p+'mlp.fc1.bias'])),w[p+'mlp.fc2.weight'],w[p+'mlp.fc2.bias'])

def encode(tokens,w):
    x=w['embeddings.word_embeddings.weight'][tokens]
    for i in range(2):x=block(x,w,i)
    return norm(x,w['ln_f.weight'],w['ln_f.bias'])

def pair_features(ref,alt,pos,w):
    B=ref.shape[0];h=encode(mx.concatenate([ref,alt],axis=0),w);a,b=h[:B],h[B:]
    diff=mx.abs(a-b);local=mx.abs(a[mx.arange(B),pos]-b[mx.arange(B),pos])
    return mx.concatenate([mx.mean((a+b)/2,axis=1),mx.mean(diff,axis=1),mx.max(diff,axis=1),local],axis=1)

def strand_features(ref,alt,pos,w):
    return .5*(pair_features(ref,alt,pos,w)+pair_features(17-ref[:,::-1],17-alt[:,::-1],ref.shape[1]-1-pos,w))

def head(features,p):return mx.tanh(features@p['head.w1']+p['head.b1'])@p['head.w2']+p['head.b2']

def bce(z,y,positive_weight):return mx.mean((1-y)*mx.logaddexp(0,z)+positive_weight*y*mx.logaddexp(0,-z))

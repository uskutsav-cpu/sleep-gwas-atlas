"""Research NumPy/Autograd port of the pinned public tiny HyenaDNA encoder.

Formula source: LongSafari/hyenadna-tiny-1k-seqlen-hf, commit e8c1effa.
Autograd supplies reverse-mode derivatives through all trained operations.
This is not an official LongSafari runtime; see port_check.json.
"""
from pathlib import Path
import sys,json,struct
sys.path.insert(0,str(Path(__file__).resolve().parent/'runtime'))
import numpy as onp
import autograd.numpy as np

WORK=Path(__file__).resolve().parent

def read_weights():
    path=WORK/'hyenadna/model.safetensors'
    with path.open('rb') as f:
        size=struct.unpack('<Q',f.read(8))[0];header=json.loads(f.read(size));raw=f.read()
    result={}
    for key,item in header.items():
        if key=='__metadata__':continue
        assert item['dtype']=='F32'
        start,end=item['data_offsets'];a=onp.frombuffer(raw[start:end],dtype='<f4').reshape(item['shape']).copy()
        result[key.removeprefix('hyena.backbone.')]=a
    return result

def linear(x,w,b=None):
    y=np.dot(x,np.transpose(w))
    return y if b is None else y+b

def norm(x,w,b):
    m=np.mean(x,axis=-1,keepdims=True)
    v=np.mean((x-m)**2,axis=-1,keepdims=True)
    return (x-m)/np.sqrt(v+1e-5)*w+b

def gelu(x):return .5*x*(1+np.tanh(.7978845608028654*(x+.044715*x**3)))

def kernel(w,p,L):
    z=w[p+'filter_fn.pos_emb.z'][:,:L,:]
    t=w[p+'filter_fn.pos_emb.t'][:,:L,:]
    f=p+'filter_fn.implicit_filter.'
    freq=w[f+'1.freq']
    h=np.sin(linear(z,w[f+'0.weight'],w[f+'0.bias'])*freq)
    for i in [2,4]:h=np.sin(linear(h,w[f+f'{i}.weight'],w[f+f'{i}.bias'])*freq)
    h=linear(h,w[f+'6.weight'])
    return h*(np.exp(-t*np.abs(w[p+'filter_fn.modulation.deltas']))+.05)

def fftconv(u,k,d):
    # Equivalent to the official /fft_size and irfft(norm='forward') pair.
    L=u.shape[-1]
    y=np.fft.irfft(np.fft.rfft(u,2*L)*np.fft.rfft(k,2*L),2*L)[...,:L]
    return y+u*d[:,None]

def direct_conv(u,k,d):
    L=u.shape[-1]
    # Independent convolution-theorem check for the long-filter implementation.
    return onp.stack([onp.stack([onp.convolve(a,b)[:L] for a,b in zip(sample,k)]) for sample in u])+u*d[:,None]

def block(x,w,i):
    p=f'layers.{i}.';m=p+'mixer.'
    h=norm(x,w[p+'norm1.weight'],w[p+'norm1.bias'])
    h=linear(h,w[m+'in_proj.weight'],w[m+'in_proj.bias'])
    B,L,D=h.shape
    pad=np.concatenate([np.zeros((B,2,D),dtype=h.dtype),h],axis=1)
    filt=w[m+'short_filter.weight'][:,0,:]
    uc=sum(pad[:,j:j+L,:]*filt[:,j] for j in range(3))+w[m+'short_filter.bias']
    x0,x1,v=np.split(uc,3,axis=-1)
    k=kernel(w,m,L)[0].T
    v=np.transpose(v*x1,(0,2,1))
    y=np.transpose(fftconv(v,k,w[m+'filter_fn.bias']),(0,2,1))*x0
    x=x+linear(y,w[m+'out_proj.weight'],w[m+'out_proj.bias'])
    h=norm(x,w[p+'norm2.weight'],w[p+'norm2.bias'])
    h=gelu(linear(h,w[p+'mlp.fc1.weight'],w[p+'mlp.fc1.bias']))
    return x+linear(h,w[p+'mlp.fc2.weight'],w[p+'mlp.fc2.bias'])

def encode(tokens,w):
    x=w['embeddings.word_embeddings.weight'][tokens]
    for i in range(2):x=block(x,w,i)
    return norm(x,w['ln_f.weight'],w['ln_f.bias'])

def pair_features(ref,alt,pos,w):
    batch=ref.shape[0]
    h=encode(onp.concatenate([ref,alt],axis=0),w)
    a,b=h[:batch],h[batch:]
    diff=np.abs(a-b)
    local=np.abs(a[onp.arange(batch),pos]-b[onp.arange(batch),pos])
    return np.concatenate([np.mean((a+b)/2,axis=1),np.mean(diff,axis=1),np.max(diff,axis=1),local],axis=1)

def strand_features(ref,alt,pos,w):
    L=ref.shape[1]
    return .5*(pair_features(ref,alt,pos,w)+pair_features(17-ref[:,::-1],17-alt[:,::-1],L-1-pos,w))

def head(features,p):
    h=np.tanh(np.dot(features,p['head.w1'])+p['head.b1'])
    return np.dot(h,p['head.w2'])+p['head.b2']

def init_head(dim,prevalence,seed=42):
    r=onp.random.default_rng(seed)
    return {'head.w1':(r.normal(size=(dim,32))*.02).astype('float32'),'head.b1':onp.zeros(32,dtype='float32'),'head.w2':(r.normal(size=32)*.02).astype('float32'),'head.b2':onp.array(onp.log(prevalence/(1-prevalence)),dtype='float32')}

def sigmoid(z):return 1/(1+onp.exp(-onp.clip(z,-30,30)))

def bce(logits,y,positive_weight):
    # Stable softplus without data-dependent mutation; this is an uncalibrated
    # detected-call ranking loss when positive_weight differs from one.
    return np.mean((1-y)*np.logaddexp(0,logits)+positive_weight*y*np.logaddexp(0,-logits))

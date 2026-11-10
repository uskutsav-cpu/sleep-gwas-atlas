"""Variant-centered51bp features plus native context, without measured outcomes."""
import sys,collections
from common import *
sys.path.insert(0,str(OLD));sys.path.insert(0,str(OLD/'runtime'))
import model,model_mlx,mlx.core as mx
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v3'))
import motif_model
LUT=np.zeros(256,dtype='uint8');LUT[[ord(c) for c in 'ACGT']]=[7,8,9,10]
SLICES={'native_context':[0,128],'variant51_hyena':[128,512],'pattern51':[512,680],'motif51':[680,1145]}
def tokenize(s):return LUT[np.frombuffer(s.encode(),dtype='uint8')]
def patterns(ref,alt):
    def counts(a):
        a=a.astype('int32')-7;parts=[]
        for k in [1,2,3]:
            ids=sum(a[:,o:a.shape[1]-k+1+o]*(4**(k-o-1)) for o in range(k));out=np.zeros((len(a),4**k),dtype='float32');np.add.at(out,(np.arange(len(a))[:,None],ids),1);parts.append(out)
        return np.concatenate(parts,1)
    def one(a,b):
        x,y=counts(a),counts(b);den=np.concatenate([np.full(4**k,a.shape[1]-k+1) for k in [1,2,3]])
        return np.c_[.5*(x+y)/den,np.abs(x-y)]
    return (.5*(one(ref,alt)+one(17-ref[:,::-1],17-alt[:,::-1]))).astype('float32')
class Encoder:
    def __init__(self):
        receipt=json.loads((V6OUT/'fixed_features_receipt.json').read_text());assert sha(OLD/'hyenadna/model.safetensors')==receipt['source_base_sha256']
        selected=json.loads((V6OUT/'selected_motifs.json').read_text());assert sha(ROOT/'work/sleep_frailty_ml_v3/motifs/H13CORE_pwms.txt')==selected['source_sha256']
        flat={};key=None
        for line in (ROOT/'work/sleep_frailty_ml_v3/motifs/H13CORE_pwms.txt').read_text().splitlines():
            if line.startswith('>'):key=line[1:];flat[key]=[]
            elif line.strip():flat[key].append([float(x) for x in line.split()])
        groups=collections.defaultdict(list)
        for r in selected['motifs']:groups[r['length']].append(np.array(flat[r['motif']],dtype='float32'))
        self.numpy_motifs={w:np.stack(g) for w,g in sorted(groups.items())};self.motifs={w:mx.array(p) for w,p in self.numpy_motifs.items()};self.labels=selected['motifs']
        assert sum(len(g) for g in self.numpy_motifs.values())==155
        self.numpy_weights=model.read_weights();self.weights={k:mx.array(v) for k,v in self.numpy_weights.items()}
        mx.metal.set_memory_limit(1024*1024**2);mx.metal.set_cache_limit(64*1024**2)
        self.encode=mx.compile(lambda a,b,p:model_mlx.strand_features(a,b,p,self.weights));self.scan=mx.compile(lambda a,b,p:motif_model.scan(a,b,p,self.motifs))
    def features(self,ref,alt,pos,context=None,batch_size=64):
        n=len(ref);result=np.empty((n,1145),dtype='float32');a51=[];b51=[];p51=[];ca=[];cb=[];cp=[];bylength=collections.defaultdict(list)
        for i,(a,b,p) in enumerate(zip(ref,alt,pos)):
            assert len(a)==len(b) and 51<=len(a)<=1024 and set(a+b)<=set('ACGT') and sum(x!=y for x,y in zip(a,b))==1 and a[p]!=b[p]
            start=max(0,min(p-25,len(a)-51));a51.append(tokenize(a[start:start+51]));b51.append(tokenize(b[start:start+51]));p51.append(p-start)
            if context is None:
                width=min(151,len(a));start=max(0,min(p-75,len(a)-width));ca.append(tokenize(a[start:start+width]));cb.append(tokenize(b[start:start+width]));cp.append(p-start);bylength[width].append(i)
        a51=np.stack(a51);b51=np.stack(b51);p51=np.array(p51,dtype='int32')
        for start in range(0,n,batch_size):
            end=min(start+batch_size,n);a=mx.array(a51[start:end].astype('int32'));b=mx.array(b51[start:end].astype('int32'));p=mx.array(p51[start:end]);result[start:end,128:512]=np.array(self.encode(a,b,p))[:,128:];result[start:end,512:680]=patterns(a51[start:end],b51[start:end]);result[start:end,680:]=np.array(self.scan(a,b,p))
        if context is not None:
            assert context.shape==(n,128);result[:,:128]=context
        else:
            for width,ids in bylength.items():
                for start in range(0,len(ids),batch_size):
                    ix=ids[start:start+batch_size];a=mx.array(np.stack([ca[i] for i in ix]).astype('int32'));b=mx.array(np.stack([cb[i] for i in ix]).astype('int32'));p=mx.array(np.array([cp[i] for i in ix],dtype='int32'));result[ix,:128]=np.array(self.encode(a,b,p))[:,:128]
        assert np.isfinite(result).all();return result
def branches():
    context=np.r_[np.arange(128),np.arange(512,596),np.arange(680,835),np.arange(990,1145)].astype('int32');variant=np.r_[np.arange(128,512),np.arange(596,680),np.arange(835,990)].astype('int32')
    assert len(context)==522 and len(variant)==623 and len(set(context)|set(variant))==1145
    return context,variant

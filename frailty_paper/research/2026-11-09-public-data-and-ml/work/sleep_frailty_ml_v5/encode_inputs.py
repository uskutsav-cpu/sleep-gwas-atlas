"""Raw paired-sequence features exactly matching the fixed V5 training bank."""
import sys,collections
from common import *
from prepare_features import selected_motifs,pattern,model,model_mlx,mx,motif_model
LUT=np.zeros(256,dtype='uint8');LUT[[ord(c) for c in 'ACGT']]=[7,8,9,10]
def tok(s):return LUT[np.frombuffer(s.encode(),dtype='uint8')]
def features(ref,alt,pos,batch_size=32):
    receipt=json.loads((OUT/'fixed_features_receipt.json').read_text());assert receipt['source_base_sha256']==sha(OLD/'hyenadna/model.safetensors');groups,labels=selected_motifs();assert labels==json.loads((OUT/'selected_motifs.json').read_text())['motifs'];D=receipt['feature_count'];weights={k:mx.array(v) for k,v in model.read_weights().items()};groups={k:mx.array(v) for k,v in groups.items()};mx.metal.set_memory_limit(1536*1024**2);mx.metal.set_cache_limit(64*1024**2);encode=mx.compile(lambda a,b,p:model_mlx.strand_features(a,b,p,weights));scan=mx.compile(lambda a,b,p:motif_model.scan(a,b,p,groups));result=np.empty((len(ref),D),dtype='float32');lr=[];la=[];lp=[];bylen=collections.defaultdict(list)
    for i,(a,b,p) in enumerate(zip(ref,alt,pos)):
        assert len(a)==len(b) and 151<=len(a)<=1024 and set(a+b)<=set('ACGT') and 0<=p<len(a) and a[p]!=b[p] and sum(c!=d for c,d in zip(a,b))==1
        start=max(0,min(p-75,len(a)-151));lr.append(tok(a[start:start+151]));la.append(tok(b[start:start+151]));lp.append(p-start);bylen[len(a)].append(i)
    lr=np.stack(lr);la=np.stack(la);lp=np.array(lp,dtype='int32')
    for start in range(0,len(ref),batch_size):
        end=min(start+batch_size,len(ref));a=mx.array(lr[start:end].astype('int32'));b=mx.array(la[start:end].astype('int32'));p=mx.array(lp[start:end]);result[start:end,:512]=np.array(encode(a,b,p));result[start:end,1024:1192]=pattern(lr[start:end],la[start:end]);result[start:end,1192:]=np.array(scan(a,b,p))
        if start%8192==0:print('Raw V5 local sequence encodings',end,'/',len(ref),flush=True)
    for length,ids in bylen.items():
        for start in range(0,len(ids),batch_size):
            ix=ids[start:start+batch_size];a=mx.array(np.stack([tok(ref[i]) for i in ix]).astype('int32'));b=mx.array(np.stack([tok(alt[i]) for i in ix]).astype('int32'));p=mx.array(np.array([pos[i] for i in ix],dtype='int32'));result[ix,512:1024]=np.array(encode(a,b,p))
    assert np.isfinite(result).all();return result

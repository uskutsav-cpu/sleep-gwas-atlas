"""Fixed DNA features and original frozen encoder, with resumable batch shards.

No activity, effect or classification labels are features. Test sequences can be
embedded by a fixed encoder, but test labels do not fit any representation.
"""
from pathlib import Path
import sys,csv,json,hashlib,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent;OUT=ROOT/'outputs/sleep_frailty_ai_v3';OLD=ROOT/'work/sleep_frailty_ml'
sys.path.insert(0,str(OLD/'runtime'));sys.path.insert(0,str(OLD));import mlx.core as mx;import model,model_mlx
mx.metal.set_memory_limit(2*1024**3);mx.metal.set_cache_limit(128*1024**2)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def kcounts(a):
    a=a.astype('int32')-7;out=[]
    for k in [1,2,3]:
        ids=sum(a[:,o:a.shape[1]-k+1+o]*(4**(k-o-1)) for o in range(k));c=np.zeros((len(a),4**k),dtype='float32');np.add.at(c,(np.arange(len(a))[:,None],ids),1);out.append(c)
    return np.concatenate(out,1)
def dna_features(a,b):
    def features(x,y):
        p,q=kcounts(x),kcounts(y);den=np.concatenate([np.full(4**k,x.shape[1]-k+1) for k in [1,2,3]])
        # Context frequencies remove the trivial 200/270 construct-length scale.
        # Local difference counts retain their natural single-SNV scale.
        return np.c_[(p+q)/2/den,np.abs(p-q)]
    return (.5*(features(a,b)+features(17-a[:,::-1],17-b[:,::-1]))).astype('float32')
weights={k:mx.array(v) for k,v in model.read_weights().items()};fn=mx.compile(lambda a,b,p:model_mlx.strand_features(a,b,p,weights));source_hash=sha(OLD/'hyenadna/model.safetensors');receipts=[]
for study in ['salomon','siraj']:
    path=WORK/(study+'_sequences_local.npz')
    with np.load(path) as z:data={k:z[k] for k in z.files}
    ref,alt,pos=data['ref'],data['alt'],data['position'];n=len(ref);features=np.empty((n,512),dtype='float32');shards=WORK/(study+'_feature_shards');shards.mkdir(exist_ok=True)
    key=hashlib.sha256((sha(path)+source_hash+sha(OLD/'model_mlx.py')).encode()).hexdigest();start=time.time()
    if study=='salomon':
        prior=json.loads((ROOT/'outputs/sleep_frailty_ai_v2/sequence_feature_audit.json').read_text());assert prior['base_weights_sha256']==source_hash
        with np.load(OLD/'prepared.npz') as z:original={k:z[k] for k in z.files}
        ix=data['v1_row_index'];assert np.array_equal(ref,original['ref'][ix]) and np.array_equal(alt,original['alt'][ix]) and np.array_equal(pos,original['position'][ix])
        all_features=np.empty((len(original['label']),512),dtype='float32');c=original['chromosome']
        for check,indices in zip(prior['feature_cache_replay'],[np.flatnonzero(c<=16),np.flatnonzero((c>=17)&(c<=18)),np.flatnonzero(c>=19)]):
            old=OLD/check['cache'];assert sha(old)==check['sha256'];all_features[indices]=np.load(old)
        features=all_features[ix];output=WORK/(study+'_features_local.npz');np.savez_compressed(output,dna_pattern=dna_features(ref,alt),hyena_original=features)
        receipts.append({'study':study,'rows':n,'input_sha256':sha(path),'features_sha256':sha(output),'base_weights_sha256':source_hash,'signature':key,'trained_on_new_data':False,'labels_in_features':False,'reuse':'Original frozen feature caches, with full selected reference/alternate/position input equality and cache hashes verified; earlier sampled encoder replay passed'})
        print('Reused verified original Salomon embeddings',n,flush=True);continue
    # Each shard has a data/checkpoint/code signature; never reuse on name alone.
    for a in range(0,n,1024):
        b=min(n,a+1024);shard=shards/f'{a}_{b}.npz'
        if shard.exists():
            with np.load(shard) as z:
                valid=str(z['signature'])==key and z['features'].shape==(b-a,512)
                if valid:features[a:b]=z['features'];continue
        arrays=[]
        for j in range(a,b,32):
            k=min(b,j+32);arrays.append(np.array(fn(mx.array(ref[j:k].astype('int32')),mx.array(alt[j:k].astype('int32')),mx.array(pos[j:k]))))
        f=np.concatenate(arrays);assert np.isfinite(f).all();features[a:b]=f;np.savez_compressed(shard,signature=np.array(key),features=f)
        print('Original encoder',study,b,'/',n,'seconds',round(time.time()-start,1),flush=True)
    output=WORK/(study+'_features_local.npz');np.savez_compressed(output,dna_pattern=dna_features(ref,alt),hyena_original=features)
    receipts.append({'study':study,'rows':n,'input_sha256':sha(path),'features_sha256':sha(output),'base_weights_sha256':source_hash,'signature':key,'trained_on_new_data':False,'labels_in_features':False})
    mx.metal.clear_cache()
(OUT/'feature_receipt.json').write_text(json.dumps({'datasets':receipts,'scope':'Only paired sequence features and fixed original pretrained representations. Activity/skew/AVI targets not used here.','test_labels_used_for_fitting':False},indent=2)+'\n');print('FEATURE EXTRACTION COMPLETE',flush=True)

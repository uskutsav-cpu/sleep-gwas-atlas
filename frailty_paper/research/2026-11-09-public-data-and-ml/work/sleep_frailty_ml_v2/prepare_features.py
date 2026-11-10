"""Build sequence-only features; AlphaGenome joins are a separate exact-key step."""
from pathlib import Path
import sys,json,hashlib
import numpy as np

ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent;OLD=ROOT/'work/sleep_frailty_ml';OUT=ROOT/'outputs/sleep_frailty_ai_v2'
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def kmer_features(ref,alt):
    def counts(a):
        a=a.astype('int32')-7;out=[]
        for k in [1,2,3]:
            ids=sum(a[:,o:a.shape[1]-k+1+o]*(4**(k-o-1)) for o in range(k))
            c=np.zeros((len(a),4**k),dtype='float32')
            np.add.at(c,(np.arange(len(a))[:,None],ids),1)
            out.append(c)
        return np.concatenate(out,axis=1)
    def features(a,b):
        x,y=counts(a),counts(b)
        return np.concatenate([(x+y)/2,np.abs(x-y)],axis=1)
    # Symmetry removes arbitrary source orientation from the comparison.
    return .5*(features(ref,alt)+features(17-ref[:,::-1],17-alt[:,::-1]))

if __name__=='__main__':
    sys.path.insert(0,str(OLD/'runtime'));sys.path.insert(0,str(OLD))
    import mlx.core as mx
    import model,model_mlx
    data=np.load(OLD/'prepared.npz');chrom=data['chromosome'];parts=[];indices=[]
    groups=[('pretrained_train_features.npy',np.flatnonzero(chrom<=16)),('pretrained_validation_features.npy',np.flatnonzero((chrom>=17)&(chrom<=18))),('pretrained_final_test_features.npy',np.flatnonzero(chrom>=19))]
    w={k:mx.array(v) for k,v in model.read_weights().items()};fn=mx.compile(lambda a,b,p:model_mlx.strand_features(a,b,p,w))
    rng=np.random.default_rng(123);checks=[];emb=np.empty((len(chrom),512),dtype='float32')
    for name,ix in groups:
        f=np.load(OLD/name);assert f.shape==(len(ix),512) and np.isfinite(f).all()
        jj=rng.choice(len(ix),min(6,len(ix)),replace=False);kk=ix[jj]
        actual=np.array(fn(mx.array(data['ref'][kk].astype('int32')),mx.array(data['alt'][kk].astype('int32')),mx.array(data['position'][kk])))
        error=float(np.max(np.abs(actual-f[jj])));assert error<1e-3,(name,error)
        checks.append({'cache':name,'sha256':sha(OLD/name),'sampled_rows':kk.tolist(),'maximum_absolute_replay_error':error,'qualification':'Sampled replay, not full row-by-row regeneration'})
        emb[ix]=f
    km=kmer_features(data['ref'],data['alt']);assert km.shape==(len(chrom),168) and np.isfinite(km).all()
    path=WORK/'sequence_features.npz';np.savez_compressed(path,kmer=km,hyena_original=emb)
    receipt={'prepared_sha256':sha(OLD/'prepared.npz'),'manifest_sha256':sha(ROOT/'outputs/sleep_frailty_ai_v1/split_manifest.csv'),'base_weights_sha256':sha(OLD/'hyenadna/model.safetensors'),'features_sha256':sha(path),'feature_definition':'168 strand-averaged 1-3mer mean/context and absolute allele-difference counts; original pretrained 512-dimensional paired features','rows':len(chrom),'feature_cache_replay':checks,'test_note':'19-22 are a previously observed v1 benchmark, not a new independent holdout'}
    (OUT/'sequence_feature_audit.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2),flush=True)

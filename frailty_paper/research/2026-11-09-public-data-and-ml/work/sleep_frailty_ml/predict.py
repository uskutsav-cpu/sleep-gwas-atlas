"""Score measured-assay-like DNA pairs with a saved local research checkpoint."""
from pathlib import Path
import sys,argparse,csv,json,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parent/'runtime'))
import numpy as np

def predict(checkpoint,rows,backend='mlx'):
    z=np.load(checkpoint);is_cnn='conv0.w' in z.files;is_kmer='weight' in z.files and 'bias' in z.files
    if backend=='mlx':
        import mlx.core as mx
        import model_mlx as model
        if is_cnn:import model_cnn as cnn
        cast=mx.array
    else:
        import model
        cast=np.asarray
    w={k.removeprefix('encoder::'):cast(z[k]) for k in z.files if k.startswith('encoder::')};h={k:cast(z[k]) for k in z.files if k.startswith('head.')}
    if not is_cnn and not is_kmer:mean=cast(z['feature_mean']);sd=cast(z['feature_sd'])
    lut=np.zeros(256,dtype='int32');lut[[ord(x) for x in 'ACGT']]=[7,8,9,10]
    ref=[];alt=[];pos=[];seen=set()
    for r in rows:
        ident=r['variant'];assert ident not in seen,('duplicate identifier',ident);seen.add(ident)
        a,b=r['reference_sequence'].upper(),r['alternate_sequence'].upper()
        assert len(a)==len(b)==270 and not set(a+b)-set('ACGT'),('unusable sequence',ident)
        changed=[i for i,(x,y) in enumerate(zip(a,b)) if x!=y];assert len(changed)==1,('require one substitution',ident)
        ref.append(lut[np.frombuffer(a.encode(),dtype=np.uint8)]);alt.append(lut[np.frombuffer(b.encode(),dtype=np.uint8)]);pos.append(changed[0])
    if is_kmer:
        def fn(a,b,p):
            f=kmer_numpy(np.array(a),np.array(b));f=np.clip((f-z['mean'])/z['sd'],-20,20)
            return f@z['weight']+z['bias']
    elif is_cnn:
        cp={k:cast(z[k]) for k in z.files}
        if backend=='mlx':fn=mx.compile(lambda a,b,p:cnn.logits(a,b,cp))
        else:fn=lambda a,b,p:cnn_numpy(a,b,cp)
    elif backend=='mlx':fn=mx.compile(lambda a,b,p:model.head(mx.clip((model.strand_features(a,b,p,w)-mean)/sd,-20,20),h))
    else:fn=lambda a,b,p:model.head(np.clip((model.strand_features(a,b,p,w)-mean)/sd,-20,20),h)
    predictions=[]
    for i in range(0,len(rows),32):
        logits=np.array(fn(cast(np.array(ref[i:i+32])),cast(np.array(alt[i:i+32])),cast(pos[i:i+32])))
        predictions.extend(1/(1+np.exp(-np.clip(logits,-30,30))))
    assert np.isfinite(predictions).all()
    return np.array(predictions)

def kmer_numpy(ref,alt):
    rows=[]
    for a,b in zip(ref-7,alt-7):
        v=[]
        for arr in [a,b]:
            c=[]
            for k in [1,2,3]:
                ids=sum(arr[o:len(arr)-k+1+o].astype('int32')*(4**(k-o-1)) for o in range(k));c.extend(np.bincount(ids,minlength=4**k))
            v.append(np.array(c,dtype='float32'))
        rows.append(np.r_[(v[0]+v[1])/2,np.abs(v[0]-v[1])])
    return np.stack(rows)

def cnn_numpy(ref,alt,p):
    def encode(tokens):
        h=np.eye(4,dtype='float32')[tokens-7]
        for i,k in enumerate([15,5,5]):
            padded=np.pad(h,((0,0),(k//2,k//2),(0,0)))
            windows=np.lib.stride_tricks.sliding_window_view(padded,k,axis=1).transpose(0,1,3,2)
            h=np.maximum(np.einsum('blkc,okc->blo',windows,p[f'conv{i}.w'],optimize=True)+p[f'conv{i}.b'],0)
        return h
    def pair(a,b):
        B=len(a);h=encode(np.concatenate([a,b],axis=0));ra,rb=h[:B],h[B:];diff=np.abs(ra-rb)
        return np.concatenate([np.mean((ra+rb)/2,axis=1),diff.mean(axis=1),diff.max(axis=1)],axis=1)
    f=.5*(pair(ref,alt)+pair(17-ref[:,::-1],17-alt[:,::-1]))
    f=f/np.sqrt(np.mean(f*f,axis=-1,keepdims=True)+1e-6)
    return np.tanh(f@p['head.w1']+p['head.b1'])@p['head.w2']+p['head.b2']

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--pairs',required=True);p.add_argument('--checkpoint',required=True);p.add_argument('--output',required=True);p.add_argument('--backend',choices=['mlx','numpy'],default='mlx');a=p.parse_args()
    out=Path(a.output)
    if out.exists():raise FileExistsError('Refusing to overwrite an existing prediction table')
    with open(a.pairs,newline='') as f:rows=list(csv.DictReader(f))
    scores=predict(a.checkpoint,rows,a.backend);order=np.argsort(-scores);ranks=np.empty(len(scores),int);ranks[order]=np.arange(1,len(scores)+1)
    with out.open('w',newline='') as f:
        fields=['variant','chromosome','position_hg38','genome_reference_allele','tested_alternate_allele','uncalibrated_neuronal_assay_call_score','rank_within_input','interpretation'];writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for r,s,rank in zip(rows,scores,ranks):writer.writerow({**{k:r.get(k,'') for k in fields[:5]},fields[5]:float(s),fields[6]:int(rank),fields[7]:'Exploratory reporter-assay ranking; not a disease risk, causal effect, or validated result for this candidate'})
    out.with_suffix('.provenance.json').write_text(json.dumps({'checkpoint_sha256':hashlib.sha256(Path(a.checkpoint).read_bytes()).hexdigest(),'pairs_sha256':hashlib.sha256(Path(a.pairs).read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'rows':len(rows),'backend':a.backend,'scope':'Uncalibrated neuronal reporter-call score; not disease risk'},indent=2)+'\n')
    print('Wrote',len(rows),'research-only pair predictions to',out)

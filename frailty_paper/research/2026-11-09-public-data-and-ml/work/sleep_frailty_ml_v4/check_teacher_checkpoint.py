"""Independent NumPy versus MLX forward replay of the adapted encoder."""
import sys
from common import *
sys.path.insert(0,str(OLD/'runtime'));sys.path.insert(0,str(OLD))
import model,model_mlx,mlx.core as mx

def main():
    status=json.loads((OUT/'teacher_training_status.json').read_text());path=Path(status['checkpoint'])
    assert sha(path)==status['checkpoint_sha256']
    with np.load(path) as z:p={k:z[k] for k in z.files}
    w={k.removeprefix('encoder::'):v for k,v in p.items() if k.startswith('encoder::')};wm={k:mx.array(v) for k,v in w.items()}
    with np.load(DATA/'teacher_data_270.npz') as z:
        ids=np.flatnonzero(z['split']==1)[[0,94,281,1100]]
        ref=z['ref'][ids].astype('int32');alt=z['alt'][ids].astype('int32');pos=z['position'][ids]
    cpu=np.asarray(model.strand_features(ref,alt,pos,w));gpu=np.array(model_mlx.strand_features(mx.array(ref),mx.array(alt),mx.array(pos),wm))
    ferr=float(np.max(np.abs(cpu-gpu)));assert ferr<5e-5,ferr
    def head(f):
        x=np.clip((f-p['feature_mean'])/p['feature_sd'],-20,20)
        h=np.tanh(x@p['head.w1']+p['head.b1'])
        return np.c_[h@p['teacher.w']+p['teacher.b'],h@p['high.w']+p['high.b']]
    # tanh and clipping are 1-Lipschitz. Bound amplification by the actual
    # per-coordinate whitening and absolute matrix weights, plus roundoff.
    bound=np.max((np.abs(cpu-gpu)/p['feature_sd'])@np.abs(p['head.w1'])@np.abs(np.c_[p['teacher.w'],p['high.w']]))+1e-5
    perr=float(np.max(np.abs(head(cpu)-head(gpu))));assert perr<=bound,(perr,bound)
    result={'checkpoint_sha256':sha(path),'examples':len(ids),'split':'Novel teacher-validation examples, not held-out tests','encoder_max_absolute_difference':ferr,'teacher_head_max_absolute_difference':perr,'teacher_head_actual_whitening_error_bound':float(bound),'passed':True,'qualification':'Independent arithmetic replay of this local research port; not parity against the original official HyenaDNA executable.'}
    write_json(OUT/'teacher_checkpoint_replay.json',result);print(json.dumps(result))
if __name__=='__main__':main()

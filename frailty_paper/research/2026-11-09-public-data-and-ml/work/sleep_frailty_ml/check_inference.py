from pathlib import Path
import sys,csv,json,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parent/'runtime'))
import numpy as np
import mlx.core as mx
from predict import predict
import model
mx.metal.set_memory_limit(2*1024**3);mx.metal.set_cache_limit(128*1024**2)

WORK=Path(__file__).resolve().parent;ROOT=WORK.parents[1];OUT=ROOT/'outputs/sleep_frailty_ai_v1'
with (WORK/'project_sequence_pairs.csv').open() as f:rows=list(csv.DictReader(f))
examples=rows[:3];rc=str.maketrans('ACGT','TGCA');checks={}
for name in ['frozen_hyenadna.npz','finetuned_hyenadna.npz','cnn.npz']:
    file=OUT/'models'/name
    gpu=predict(file,examples,'mlx');cpu=predict(file,examples,'numpy')
    error=float(np.max(np.abs(gpu-cpu)));assert error<3e-4,(name,error)
    swapped=[{**r,'reference_sequence':r['alternate_sequence'],'alternate_sequence':r['reference_sequence']} for r in examples]
    reverse=[{**r,'reference_sequence':r['reference_sequence'].translate(rc)[::-1],'alternate_sequence':r['alternate_sequence'].translate(rc)[::-1]} for r in examples]
    swap_error=float(np.max(np.abs(gpu-predict(file,swapped,'mlx'))));rc_error=float(np.max(np.abs(gpu-predict(file,reverse,'mlx'))));assert max(swap_error,rc_error)<1e-6
    rejected=[]
    for label,bad in [('duplicate',examples+[examples[0]]),('ambiguous_base',[{**examples[0],'reference_sequence':'N'+examples[0]['reference_sequence'][1:]}]),('unchanged_pair',[{**examples[0],'alternate_sequence':examples[0]['reference_sequence']}]),('wrong_length',[{**examples[0],'reference_sequence':examples[0]['reference_sequence'][:-1]}])]:
        try:predict(file,bad,'numpy')
        except (AssertionError,ValueError):rejected.append(label)
        else:raise ValueError(('Invalid input was accepted',name,label))
    checks[name]={'cpu_gpu_prediction_max_error':error,'allele_swap_score_max_error':swap_error,'reverse_complement_score_max_error':rc_error,'invalid_inputs_rejected':rejected,'checkpoint_sha256':hashlib.sha256(file.read_bytes()).hexdigest()}

original=model.read_weights();fine=np.load(OUT/'models/finetuned_hyenadna.npz');frozen=np.load(OUT/'models/frozen_hyenadna.npz')
fixed=[k for k in original if k.startswith('layers.0.') or k.startswith('embeddings.')]
assert all(np.array_equal(fine['encoder::'+k],original[k]) for k in fixed)
assert all(np.array_equal(frozen['encoder::'+k],v) for k,v in original.items() if k!='lm_head.weight')
changed=[k for k in original if 'encoder::'+k in fine and not np.array_equal(fine['encoder::'+k],original[k])]
assert len(changed)>5 and all(k.startswith('layers.1.') or k.startswith('ln_f.') for k in changed)
receipt={'models':checks,'native_pretrained_arrays_changed':changed,'first_block_and_nucleotide_embedding_preserved':True,'frozen_encoder_matches_original':True,'scope':'Numerical reload, invariance, input rejection and parameter-update checks; not biological or clinical validation'}
# Independently reload selected checkpoints and replay their validation metric.
# This catches stale captured parameters in compiled validation functions.
data=np.load(WORK/'prepared.npz');vi=np.flatnonzero((data['chromosome']>=17)&(data['chromosome']<=18));alphabet=np.array(list('ACGT'))
validation_rows=[{'variant':str(int(i)),'reference_sequence':''.join(alphabet[data['ref'][i]-7]),'alternate_sequence':''.join(alphabet[data['alt'][i]-7])} for i in vi]
def average_precision(y,p):
    ix=np.argsort(-p,kind='stable');s=p[ix];a=y[ix];ends=np.r_[np.flatnonzero(np.diff(s)!=0),len(s)-1];tp=np.cumsum(a)[ends]
    return float(np.sum(np.diff(np.r_[0,tp/y.sum()])*tp/(ends+1)))
history=json.loads((OUT/'training_history.json').read_text());cnn=json.loads((OUT/'cnn_selection.json').read_text())
expected={'frozen_hyenadna.npz':max(x['validation_average_precision'] for x in history if x['phase']=='frozen_head'),'finetuned_hyenadna.npz':max(x['validation_average_precision'] for x in history if x['phase']=='finetune'),'cnn.npz':cnn['validation_average_precision']}
replay={}
for name,score in expected.items():
    values=predict(OUT/'models'/name,validation_rows,'mlx');actual=average_precision(data['label'][vi],values);difference=abs(actual-score)
    replay[name]={'recorded_validation_average_precision':score,'reloaded_validation_average_precision':actual,'difference':difference}
    assert difference<1e-4,(name,score,actual)
receipt['validation_checkpoint_replay']=replay
(OUT/'inference_check.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

"""Export learned project-variant representations; no interaction/causality claim."""
from pathlib import Path
import sys,csv,json,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parent/'runtime'))
import numpy as np
import mlx.core as mx
import model_mlx as model
WORK=Path(__file__).resolve().parent;ROOT=WORK.parents[1];OUT=ROOT/'outputs/sleep_frailty_ai_v1'
with (WORK/'project_sequence_pairs.csv').open() as f:rows=list(csv.DictReader(f))
checkpoint=OUT/'models/finetuned_hyenadna.npz';z=np.load(checkpoint);w={k.removeprefix('encoder::'):mx.array(z[k]) for k in z.files if k.startswith('encoder::')};h={k:mx.array(z[k]) for k in z.files if k.startswith('head.')}
lut=np.zeros(256,dtype='int32');lut[[ord(c) for c in 'ACGT']]=[7,8,9,10];ref=[];alt=[];positions=[]
for r in rows:
    a,b=r['reference_sequence'],r['alternate_sequence'];changes=[i for i,(u,v) in enumerate(zip(a,b)) if u!=v];assert len(a)==len(b)==270 and len(changes)==1
    ref.append(lut[np.frombuffer(a.encode(),dtype=np.uint8)]);alt.append(lut[np.frombuffer(b.encode(),dtype=np.uint8)]);positions.append(changes[0])
features=model.strand_features(mx.array(np.array(ref)),mx.array(np.array(alt)),mx.array(positions),w)
unclipped=(features-mx.array(z['feature_mean']))/mx.array(z['feature_sd'])
scaled=mx.clip(unclipped,-20,20)
embedding=mx.tanh(scaled@h['head.w1']+h['head.b1'])
np.savez_compressed(OUT/'project_candidate_embeddings.npz',variant_ids=np.array([r['variant'] for r in rows]),pair_features_512=np.array(features),standardized_pair_features_512=np.array(scaled),learned_head_embedding_32=np.array(embedding))
raw_z=np.array(unclipped)
with (OUT/'project_embedding_diagnostics.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=['variant','max_absolute_training_standardized_feature','fraction_features_clipped_at_20']);writer.writeheader()
    for r,zz in zip(rows,raw_z):writer.writerow({'variant':r['variant'],'max_absolute_training_standardized_feature':float(np.abs(zz).max()),'fraction_features_clipped_at_20':float(np.mean(np.abs(zz)>20))})
receipt={'variants':len(rows),'representations':{'pair_features':512,'learned_task_embedding':32},'checkpoint_sha256':hashlib.sha256(checkpoint.read_bytes()).hexdigest(),'max_absolute_training_standardized_feature':float(np.abs(raw_z).max()),'training_task':'Detected allelic effects in selected human excitatory-neuron reporter constructs','caveats':'No exact project variants were tested in the source assay. Embedding proximity is not a gene interaction, shared disease mechanism, or individual risk. These are exploratory outputs from the fine-tuned model, regardless of which model is recommended by validation. Standardized-feature/clipping diagnostics flag numerical extrapolation only; they do not establish biological in-distribution status.'}
(OUT/'project_embedding_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

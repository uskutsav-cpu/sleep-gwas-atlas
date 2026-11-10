"""Embed verified external DNA with the original frozen encoder; no fitting."""
from pathlib import Path
import sys,csv,json,hashlib
import numpy as np
from prepare_features import kmer_features
ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent;OLD=ROOT/'work/sleep_frailty_ml';OUT=ROOT/'outputs/sleep_frailty_ai_v2'
sys.path.insert(0,str(OLD/'runtime'));sys.path.insert(0,str(OLD));import mlx.core as mx;import model,model_mlx
audit=json.loads((OUT/'external_mapping_audit.json').read_text());assert not audit['retrieval_errors']
pairs=WORK/'external_pairs_local.csv';rows=list(csv.DictReader(pairs.open()));lut=np.zeros(256,dtype='int32');lut[[ord(x) for x in 'ACGT']]=[7,8,9,10]
ref=np.stack([lut[np.frombuffer(r['reference_sequence'].encode(),dtype='uint8')] for r in rows]);alt=np.stack([lut[np.frombuffer(r['alternate_sequence'].encode(),dtype='uint8')] for r in rows]);pos=np.full(len(rows),135,dtype='int32')
assert ref.shape==alt.shape==(len(rows),270) and np.all((ref!=alt).sum(1)==1)
w={k:mx.array(v) for k,v in model.read_weights().items()};fn=mx.compile(lambda a,b,p:model_mlx.strand_features(a,b,p,w));features=[]
mx.metal.set_memory_limit(2*1024**3);mx.metal.set_cache_limit(128*1024**2)
for start in range(0,len(rows),32):
    features.append(np.array(fn(mx.array(ref[start:start+32]),mx.array(alt[start:start+32]),mx.array(pos[start:start+32]))))
    if start%256==0:print('External embeddings',start,'/',len(rows),flush=True)
out=WORK/'external_features_local.npz';np.savez_compressed(out,ref=ref,alt=alt,position=pos,kmer=kmer_features(ref,alt),hyena_original=np.concatenate(features),label=np.array([int(r['label']) for r in rows],dtype='float32'))
(OUT/'external_feature_receipt.json').write_text(json.dumps({'rows':len(rows),'pairs_sha256':hashlib.sha256(pairs.read_bytes()).hexdigest(),'features_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'base_weights_sha256':hashlib.sha256((OLD/'hyenadna/model.safetensors').read_bytes()).hexdigest(),'fitting_on_external_data':False,'scope':'Original frozen DNA representations, labels reserved for external evaluation only'},indent=2)+'\n');print('External features ready',len(rows),flush=True)

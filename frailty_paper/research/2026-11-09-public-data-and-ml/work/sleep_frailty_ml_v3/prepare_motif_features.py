"""Resumable biological sequence features, stored on the existing research SSD."""
from pathlib import Path
import argparse,json,time,hashlib,sys
import numpy as np
from motif_model import load_motifs,scan,mx
ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent;OUT=ROOT/'outputs/sleep_frailty_ai_v3'
p=argparse.ArgumentParser();p.add_argument('--store',default='/Volumes/Extreme SSD/codex-sleep-frailty-work/sleep_frailty_ml_v3');args=p.parse_args();store=Path(args.store);store.mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
check=json.loads((OUT/'motif_implementation_check.json').read_text());groups,labels=load_motifs();assert len(labels)==check['motif_count'];groups={w:mx.array(v) for w,v in groups.items()};M=len(labels);fn=mx.compile(lambda a,b,p:scan(a,b,p,groups));receipts=[]
mx.metal.set_memory_limit(2*1024**3);mx.metal.set_cache_limit(128*1024**2)
source_hash=sha(WORK/'motifs/H13CORE_pwms.txt')
for stem in ['salomon','siraj']:
    path=WORK/(stem+'_sequences_local.npz')
    with np.load(path) as z:data={k:z[k] for k in ['ref','alt','position']}
    n=len(data['ref']);signature=hashlib.sha256((sha(path)+source_hash+sha(WORK/'motif_model.py')).encode()).hexdigest();target=store/(stem+'_motif_features.npy');progress=store/(stem+'_motif_progress.json')
    saved=json.loads(progress.read_text()) if progress.exists() else None
    if saved and saved.get('signature')==signature and target.exists():
        out=np.lib.format.open_memmap(target,mode='r+');assert out.shape==(n,3*M) and out.dtype==np.float16;completed=set(saved['completed_shards'])
    else:
        out=np.lib.format.open_memmap(target,mode='w+',dtype='float16',shape=(n,3*M));out[:]=np.nan;out.flush();completed=set()
    started=time.time();expected=list(range(0,n,1024))
    for start in expected:
        if start in completed:continue
        stop=min(n,start+1024)
        for a in range(start,stop,32):
            b=min(stop,a+32);f=np.array(fn(mx.array(data['ref'][a:b].astype('int32')),mx.array(data['alt'][a:b].astype('int32')),mx.array(data['position'][a:b])));assert f.shape==(b-a,3*M) and np.isfinite(f).all();out[a:b]=f.astype('float16')
        out.flush();completed.add(start);record={'signature':signature,'completed_shards':sorted(completed),'rows':n,'features':3*M,'dtype':'float16','complete':len(completed)==len(expected)};tmp=progress.with_suffix('.tmp');tmp.write_text(json.dumps(record,indent=2)+'\n');tmp.replace(progress)
        print('Motif encoder',stem,stop,'/',n,'seconds',round(time.time()-started,1),flush=True)
    assert completed==set(expected)
    for a in expected:assert np.isfinite(out[a:min(n,a+1024)]).all()
    del out
    receipts.append({'study':stem,'path':str(target),'rows':n,'feature_count':3*M,'dtype':'float16','input_sha256':sha(path),'features_sha256':sha(target),'signature':signature,'fixed_encoder':True,'test_labels_used_for_fitting':False})
    mx.metal.clear_cache()
(OUT/'motif_feature_receipt.json').write_text(json.dumps({'datasets':receipts,'motif_count':M,'feature_order':[r['motif'] for r in labels],'blocks':['mean global affinity','absolute local change','mean local affinity'],'source_pwm_sha256':source_hash,'storage':'Existing research SSD; model weights/code stay in workspace. No full feature matrix needs to be packaged.','precision':'Features computed float32 and stored float16; use float32 batches for training. Half-precision quantization remains an implementation qualification.'},indent=2)+'\n');print('MOTIF FEATURE EXTRACTION COMPLETE',flush=True)

"""Small biological encoder made from verified HOCOMOCO DNA-binding models.

Scores are sequence compatibility, not confirmed TF binding, gene regulation,
causality or a disease probability. No AlphaGenome-derived motif weights.
"""
from pathlib import Path
import sys,json,collections
import numpy as np
OLD=Path(__file__).resolve().parents[1]/'sleep_frailty_ml';WORK=Path(__file__).resolve().parent
sys.path.insert(0,str(OLD/'runtime'));import mlx.core as mx

def load_motifs(limit=None):
    meta=[json.loads(s) for s in (WORK/'motifs/H13CORE_annotation.jsonl').read_text().splitlines()]
    flat={};name=None
    for line in (WORK/'motifs/H13CORE_pwms.txt').read_text().splitlines():
        if line.startswith('>'):name=line[1:];flat[name]=[]
        elif line.strip():flat[name].append([float(x) for x in line.split()])
    selected=sorted([r for r in meta if 'HUMAN' in r['masterlist_info']['species']],key=lambda r:(r['length'],r['name']))
    if limit is not None:selected=selected[:limit]
    groups=collections.defaultdict(list);labels=[]
    for r in selected:
        pwm=np.asarray(flat[r['name']],dtype='float32');assert pwm.shape==(r['length'],4) and np.isfinite(pwm).all()
        assert np.max(np.abs(pwm-np.asarray(r['pwm'])))<1e-6
        groups[r['length']].append(pwm);labels.append({'motif':r['name'],'TF':r['tf'],'length':r['length'],'quality':r['quality'],'human_gene':r['masterlist_info']['species']['HUMAN']['gene_symbol']})
    return {w:np.stack(v) for w,v in sorted(groups.items())},labels

def scan(ref,alt,pos,groups):
    B=ref.shape[0];tokens=mx.concatenate([ref,alt],0);position=mx.concatenate([pos,pos],0);onehot=mx.eye(4)[tokens-7];global_features=[];local_features=[]
    for width,pwm in groups.items():
        n=pwm.shape[0];both=mx.concatenate([pwm,pwm[:,::-1,::-1]],0);raw=mx.conv1d(onehot,both)
        score=mx.maximum(raw[:,:,:n],raw[:,:,n:]);lo=mx.min(pwm,axis=2).sum(axis=1);hi=mx.max(pwm,axis=2).sum(axis=1)
        score=mx.clip((score-lo)/(hi-lo),0,1);global_features.append(mx.max(score,axis=1))
        start=mx.arange(ref.shape[1]-width+1);overlap=(start[None,:]<=position[:,None]) & (position[:,None]<start[None,:]+width)
        local_features.append(mx.max(mx.where(overlap[:,:,None],score,-1e9),axis=1))
    g=mx.concatenate(global_features,1);loc=mx.concatenate(local_features,1);ga,gb=g[:B],g[B:];la,lb=loc[:B],loc[B:]
    return mx.concatenate([.5*(ga+gb),mx.abs(la-lb),.5*(la+lb)],1)

def scan_numpy(ref,alt,pos,groups):
    B=len(ref);tokens=np.concatenate([ref,alt]);position=np.r_[pos,pos];onehot=np.eye(4,dtype='float32')[tokens-7];glob=[];local=[]
    for width,pwm in groups.items():
        both=np.concatenate([pwm,pwm[:,::-1,::-1]],0);windows=np.lib.stride_tricks.sliding_window_view(onehot,width,axis=1).transpose(0,1,3,2)
        raw=np.einsum('blwc,mwc->blm',windows,both,optimize=True);n=len(pwm);score=np.maximum(raw[:,:,:n],raw[:,:,n:]);lo=pwm.min(2).sum(1);hi=pwm.max(2).sum(1);score=np.clip((score-lo)/(hi-lo),0,1);glob.append(score.max(1))
        start=np.arange(ref.shape[1]-width+1);overlap=(start[None,:]<=position[:,None]) & (position[:,None]<start[None,:]+width);local.append(np.where(overlap[:,:,None],score,-1e9).max(1))
    g=np.concatenate(glob,1);loc=np.concatenate(local,1)
    return np.c_[.5*(g[:B]+g[B:]),np.abs(loc[:B]-loc[B:]),.5*(loc[:B]+loc[B:])]

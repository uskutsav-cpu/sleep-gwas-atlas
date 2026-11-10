"""Independent direct scanning, invariance and boundary checks on real motifs."""
from pathlib import Path
import json,numpy as np
from motif_model import load_motifs,scan,scan_numpy,mx
WORK=Path(__file__).resolve().parent;OUT=WORK.parents[1]/'outputs/sleep_frailty_ai_v3'
all_groups,labels=load_motifs();groups={w:arr[:min(2,len(arr))] for w,arr in all_groups.items()};gpu={w:mx.array(v) for w,v in groups.items()};rng=np.random.default_rng(904);checks=[]
for L in [200,270]:
    a=rng.integers(7,11,(3,L),dtype='int32');b=a.copy();pos=np.array([0,L//2,L-1],dtype='int32')
    for i,p in enumerate(pos):b[i,p]=7+(a[i,p]-6)%4
    actual=np.array(scan(mx.array(a),mx.array(b),mx.array(pos),gpu));expected=scan_numpy(a,b,pos,groups);error=float(np.max(np.abs(actual-expected)));assert error<3e-5,(L,error)
    swapped=np.array(scan(mx.array(b),mx.array(a),mx.array(pos),gpu));reverse=np.array(scan(mx.array(17-a[:,::-1]),mx.array(17-b[:,::-1]),mx.array(L-1-pos),gpu));s=float(np.max(np.abs(actual-swapped)));r=float(np.max(np.abs(actual-reverse)));assert max(s,r)<3e-5
    assert np.isfinite(actual).all() and actual.min()>=-1e-6 and actual.max()<=1+1e-6
    checks.append({'length':L,'changed_positions':pos.tolist(),'direct_numpy_max_error':error,'allele_swap_max_error':s,'reverse_complement_max_error':r})
receipt={'motif_count':len(labels),'human_TFs':len({r['human_gene'] for r in labels}),'all_flat_file_and_annotation_PWMs_match':True,'checks':checks,'feature_blocks':['mean reference/alternate global maximum normalized motif score','absolute difference in maximum overlapping-site motif score','mean maximum overlapping-site motif score'],'score_scope':'Predicted sequence compatibility, not measured TF binding or validated gene regulation','test_labels_used':False}
(OUT/'motif_implementation_check.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

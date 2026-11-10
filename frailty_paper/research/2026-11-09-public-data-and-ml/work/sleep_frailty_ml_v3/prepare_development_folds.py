"""Audit exact 200-base DNA matches across three development-only folds.

No model inference and no inspection of reserved chromosomes1–3 outcomes.
Full SHA256 digests retain all valid mutation-containing windows and both alleles.
"""
from common import *
FOLDS=[list(range(4,10)),list(range(10,16)),list(range(16,23))]
source={};arrays=[];offset=0;rows=[]
for stem in ['salomon','siraj']:
    p=WORK/(stem+'_sequences_local.npz');source[stem]=sha(p)
    with np.load(p) as z:
        a={k:z[k] for k in ['ref','alt','position','chromosome','label']}
    for j in np.flatnonzero(a['chromosome']>=4):rows.append((offset+int(j),int(a['chromosome'][j]),a['ref'][j],a['alt'][j],int(a['position'][j])))
    arrays.append(a);offset+=len(a['label'])
chrom=np.concatenate([a['chromosome'] for a in arrays]);study=np.concatenate([np.full(len(a['label']),i) for i,a in enumerate(arrays)]);y=np.concatenate([a['label'] for a in arrays])
def windows(a,p):
    assert len(a)>=200 and 0<=p<len(a)
    # Original encoder tokens:7=A,8=C,9=G,10=T.
    assert np.all((a>=7)&(a<=10))
    for start in range(max(0,p-199),min(p,len(a)-200)+1):
        w=a[start:start+200];b=w.tobytes();r=(17-w[::-1]).tobytes();yield hashlib.sha256(min(b,r)).digest()
assert np.array_equal(17-np.array([7,8,9,10],dtype='uint8')[::-1],np.array([7,8,9,10],dtype='uint8'))
result=[]
for val in FOLDS:
    held=set()
    for i,c,a,b,p in rows:
        if c in val:held.update(windows(a,p));held.update(windows(b,p))
    excluded=[]
    for i,c,a,b,p in rows:
        if c not in val and (any(w in held for w in windows(a,p)) or any(w in held for w in windows(b,p))):excluded.append(i)
    tr=np.flatnonzero((chrom>=4)&~np.isin(chrom,val));tr=tr[~np.isin(tr,excluded)];va=np.flatnonzero(np.isin(chrom,val))
    r={'validation_chromosomes':val,'excluded_training_indices':excluded,'training_indices_sha256':hashlib.sha256(tr.tobytes()).hexdigest(),'validation_indices_sha256':hashlib.sha256(va.tobytes()).hexdigest(),'training_n':len(tr),'validation_n':len(va),'studies':{STUDIES[s]:{'training_n':int(np.sum(study[tr]==s)),'training_calls':int(y[tr[study[tr]==s]].sum()),'validation_n':int(np.sum(study[va]==s)),'validation_calls':int(y[va[study[va]==s]].sum())} for s in range(2)}}
    result.append(r);print('Audited development fold',val,'removed',len(excluded),'training',len(tr),'validation',len(va),flush=True)
out={'source_sha256':source,'folds':result,'window_rule':'SHA256 of canonical full200base windows containing the tested substitution, both alleles; exclude training matches to validation','test_chromosomes_not_used':[1,2,3],'scope':'Development robustness checks after initial4–5 screening; not unbiased nested CV or independent biological replication'}
(OUT/'development_fold_sequence_audit.json').write_text(json.dumps(out,indent=2)+'\n')

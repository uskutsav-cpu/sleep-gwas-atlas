#!/usr/bin/env python3
import gzip, hashlib, json
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent
with gzip.open(p/'strand_safe_ld.tsv.gz','rt') as f: R=np.loadtxt(f,delimiter='\t')
assert np.isfinite(R).all() and np.max(np.abs(R-R.T)) < 1e-7
vals, vecs=np.linalg.eigh(R)
clip=(vecs * np.maximum(vals,0)) @ vecs.T
d=np.diag(clip).copy(); assert np.all(d>0)
clip=clip/np.sqrt(np.outer(d,d))
ridge=.999*R + .001*np.eye(len(R))
out={}
for name,A in [('eigenclip',clip),('ridge_0p001',ridge)]:
    eig=np.linalg.eigvalsh(A)
    delta=A-R
    report={'min_eigenvalue':float(eig.min()),'max_asymmetry':float(np.max(abs(A-A.T))),'max_diagonal_deviation':float(np.max(abs(np.diag(A)-1))),'max_abs_entry_change':float(np.max(abs(delta))),'relative_frobenius_change':float(np.linalg.norm(delta)/np.linalg.norm(R))}
    assert report['min_eigenvalue'] > -1e-8 and report['max_asymmetry'] < 1e-8 and report['max_diagonal_deviation'] < 1e-7
    file=p/(name+'.tsv.gz')
    with gzip.open(file,'wt') as f: np.savetxt(f,A,delimiter='\t',fmt='%.12g')
    report['sha256']=hashlib.sha256(file.read_bytes()).hexdigest();out[name]=report
(p/'repair_sensitivity_receipt.json').write_text(json.dumps({'protocol':'REPAIR_SENSITIVITY_PROTOCOL.md','input_min_eigenvalue':float(vals.min()),'transforms':out,'note':'Post-hoc numerical sensitivity only; neither transformation validates study LD.'},indent=2)+'\n')
print(json.dumps(out,indent=2),flush=True)

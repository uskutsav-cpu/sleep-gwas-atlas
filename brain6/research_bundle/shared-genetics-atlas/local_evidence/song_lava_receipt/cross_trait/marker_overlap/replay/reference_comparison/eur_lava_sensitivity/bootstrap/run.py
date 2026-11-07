import json,hashlib
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent;src=p.parent/'eur.py';source=src.read_text()
# Reuse the same documented public-data retrieval and baseline computation, but
# stop before its output-writing section so existing evidence is not overwritten.
prefix=source.split("with (p/'support.tsv').open('w') as f:",1)[0]
scope={'__file__':str(src),'__name__':'reference_bootstrap_reconstruction'}
exec(compile(prefix,str(src),'exec'),scope)
X=scope['X'];selected=scope['selected'];baseline=json.loads((p.parent/'eur_result.json').read_text())
assert X.shape==(503,2738) and len(selected)==baseline['matched_markers']
assert abs(scope['omega']-baseline['fi_omega'])<1e-12 and scope['keep']==baseline['K_99']
z=np.array([float(m['STAT']) for m in selected]);n=np.array([float(m['N']) for m in selected]);N=float(n.mean())
corr=z/np.sqrt(z*z+n-2);rng=np.random.default_rng(46420261028);results=[]
for iteration in range(1,31):
    draw=rng.integers(0,503,size=503);xb=X[draw,:];xb=xb-xb.mean(axis=0);sd=xb.std(axis=0,ddof=1)
    ok=sd>1e-12;xb=xb[:,ok]/sd[ok];U,d,Vt=np.linalg.svd(xb,full_matrices=False)
    lam=d*d/502;k=int(np.flatnonzero(np.cumsum(lam)/max(np.cumsum(lam))*100>=99)[0]+1)
    delta=(Vt[:k,:]@corr[ok])/np.sqrt(lam[:k]);dtd=float(delta@delta);sigma=(1-dtd)/(N-k-1);omega=dtd/k-sigma
    results.append({'replicate':iteration,'n_snps':int(ok.sum()),'K':k,'delta_sumsq':dtd,'sigma':sigma,'fi_omega':omega,'fi_h2_uncapped':dtd-k*sigma})
    print('replicate',iteration,'K',k,'omega',omega,flush=True)
out={'seed':46420261028,'replicates':results,'baseline_fi_omega':baseline['fi_omega'],
     'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
     'note':'Resampling one small reference, not independent cohorts or calibrated confidence intervals; no genotypes saved.'}
(p/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print('BOOTSTRAP_COMPLETE',flush=True)

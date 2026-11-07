import csv,json,hashlib
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent;src=p.parent/'eur.py';source=src.read_text()
# Run only the parent's public-data loading and baseline computation; preserve
# existing files and keep the per-person genotype matrix in memory.
prefix=source.split("with (p/'support.tsv').open('w') as f:",1)[0]
scope={'__file__':str(src),'__name__':'reference_population_holdout'}
exec(compile(prefix,str(src),'exec'),scope)
X=scope['X'];selected=scope['selected'];samples=scope['samples'];baseline=json.loads((p.parent/'eur_result.json').read_text())
assert X.shape==(503,2738) and abs(scope['omega']-baseline['fi_omega'])<1e-12 and scope['keep']==baseline['K_99']
panel=list(csv.DictReader(scope['panel'].open(),delimiter='\t'));by_id={r['sample']:r['pop'] for r in panel if r['super_pop']=='EUR'}
pops=np.array([by_id[s] for s in samples]);groups=sorted(set(pops));assert groups==['CEU','FIN','GBR','IBS','TSI']
z=np.array([float(m['STAT']) for m in selected]);n=np.array([float(m['N']) for m in selected]);N=float(n.mean());corr=z/np.sqrt(z*z+n-2)
rows=[]
for excluded in groups:
 xb=X[pops!=excluded,:];nref=xb.shape[0];xb=xb-xb.mean(axis=0);sd=xb.std(axis=0,ddof=1);ok=sd>1e-12
 xb=xb[:,ok]/sd[ok];U,d,Vt=np.linalg.svd(xb,full_matrices=False);lam=d*d/(nref-1)
 k=int(np.flatnonzero(np.cumsum(lam)/max(np.cumsum(lam))*100>=99)[0]+1)
 delta=(Vt[:k,:]@corr[ok])/np.sqrt(lam[:k]);dtd=float(delta@delta);sigma=(1-dtd)/(N-k-1);omega=dtd/k-sigma
 rows.append({'excluded_population':excluded,'n_people':nref,'n_snps':int(ok.sum()),'K':k,'delta_sumsq':dtd,'sigma':sigma,'fi_omega':omega,'fi_h2_uncapped':dtd-k*sigma})
 print(excluded,'people',nref,'snps',ok.sum(),'K',k,'omega',omega,flush=True)
out={'baseline_omega':baseline['fi_omega'],'baseline_K':baseline['K_99'],'holdouts':rows,
 'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
 'note':'Overlapping reference subsets; not independent cohorts, calibrated inference, or a local-correlation result.'}
(p/'results.json').write_text(json.dumps(out,indent=2)+'\n');print('HOLDOUT_COMPLETE',flush=True)

#!/usr/bin/env python3
import csv, gzip, json
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent
old=p.parents[2]/'rs62264764_dense_qtl/multisignal_pilot'
with (p/'matched_snps.tsv').open() as f: matched=list(csv.DictReader(f,delimiter='\t'))
with (old/'variants.tsv').open() as f: original=list(csv.DictReader(f,delimiter='\t'))
lookup={v['variant']:v for v in original}
data=np.load(p/'matched_ld.npz'); names=data['variant']; R=data['ld']
assert list(names)==[x['variant'] for x in matched]
keep=np.array([set((x['ref'],x['alt'])) not in ({'A','T'},{'C','G'}) for x in matched])
names=names[keep]; rows=[x for x,k in zip(matched,keep) if k]; R=R[np.ix_(keep,keep)]
freq=np.array([float(x['gctb_a1freq']) if x['gctb_a1']==x['alt'] else 1-float(x['gctb_a1freq']) for x in rows])
oldfreq=np.array([float(lookup[x['variant']]['reference_alt_frequency']) for x in rows])
diff=np.abs(freq-oldfreq)
assert np.isfinite(R).all() and np.max(np.abs(R-R.T))<1e-6
diag=np.diag(R).copy(); assert np.all(diag>0)
R=R / np.sqrt(np.outer(diag,diag))
eig=np.linalg.eigvalsh(R.astype('float64'))
report={'matched':len(matched),'strand_ambiguous_excluded':int((~keep).sum()),'strand_safe':len(rows),'raw_diagonal_min':float(diag.min()),'raw_diagonal_max':float(diag.max()),'normalization':'R_ij divided by sqrt(R_ii * R_jj) to obtain unit diagonal','correlation_min':float(R.min()),'correlation_max':float(R.max()),'minimum_eigenvalue':float(eig.min()),'eigenvalues_below_minus_1e-6':int((eig < -1e-6).sum()),'frequency_difference_median':float(np.median(diff)),'frequency_difference_max':float(diff.max()),'frequency_difference_over_0_1':int((diff>.1).sum()),'limitations':'Reference allele frequencies compared with 503-person reference; both are external to QTL and GWAS. Numerical checks do not prove LD–Z compatibility.'}
with gzip.open(p/'strand_safe_ld.tsv.gz','wt') as f:
    np.savetxt(f,R,delimiter='\t',fmt='%.9g')
(p/'strand_safe_variants.txt').write_text('\n'.join(names)+'\n')
(p/'qc_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

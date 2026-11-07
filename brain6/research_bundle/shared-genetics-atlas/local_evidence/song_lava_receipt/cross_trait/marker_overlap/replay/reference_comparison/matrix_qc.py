import csv,json,hashlib
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent
with (p/'ukb_ld.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
ids=[r['SNP'] for r in rows];a=np.array([[float(r[s]) for s in ids] for r in rows])
assert a.shape==(60,60) and np.all(np.isfinite(a)) and np.max(abs(a-a.T))<1e-12
assert np.allclose(np.diag(a),1,atol=1e-12)
e=np.linalg.eigvalsh(a);i,j=np.triu_indices(len(ids),1); bad=[{'snp1':ids[x],'snp2':ids[y],'value':float(a[x,y])} for x,y in zip(i,j) if abs(a[x,y])>1]
out={'dimension':60,'minimum_eigenvalue':float(e[0]),'negative_eigenvalues_below_minus_1e-8':int(sum(e < -1e-8)),
     'offdiagonal_abs_over_one':bad,'ukb_matrix_sha256':hashlib.sha256((p/'ukb_ld.tsv').read_bytes()).hexdigest(),
     'note':'sampled 60-marker approximate LD matrix only; does not establish the cause of the locus-wide LAVA failure'}
(p/'matrix_qc.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

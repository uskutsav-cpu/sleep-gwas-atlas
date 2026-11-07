import json,subprocess,hashlib
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent;src=p.parent/'eur.py';source=src.read_text()
prefix=source.split("with (p/'support.tsv').open('w') as f:",1)[0]
scope={'__file__':str(src),'__name__':'reference_runtime_reconstruction'}
exec(compile(prefix,str(src),'exec'),scope)
baseline=json.loads((p.parent/'eur_result.json').read_text());X=scope['X'];selected=scope['selected']
assert X.shape==(503,2738) and scope['keep']==baseline['K_99'] and abs(scope['omega']-baseline['fi_omega'])<1e-12
with (p.parent/'support.tsv').open() as f:
 import csv
 assert [r['SNP'] for r in csv.DictReader(f,delimiter='\t')]==[r['SNP'] for r in selected]
ld=(X.T@X)/502
assert np.max(np.abs(ld-ld.T))<1e-12 and np.max(np.abs(np.diag(ld)-1))<1e-10
cmd=['/private/tmp/brain6ext/r-env-4.3.3-short/bin/Rscript',str(p/'run.R')]
r=subprocess.run(cmd,input=np.asarray(ld,dtype='<f8').tobytes(order='F'),capture_output=True,timeout=180)
assert r.returncode==0,(r.returncode,r.stdout.decode()[-2000:],r.stderr.decode()[-2000:])
assert 'LAVA_UNIVARIATE_CHECK_COMPLETE' in r.stdout.decode()
(p/'receipt.json').write_text(json.dumps({'input_ld_shape':list(ld.shape),'baseline_source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
 'R_exit_code':r.returncode,'R_stdout':r.stdout.decode(),'R_stderr':r.stderr.decode(),
 'note':'LD streamed in memory; original UKB allele frequencies retained; no bivariate test or genotype files.'},indent=2)+'\n')
print(r.stdout.decode(),flush=True)

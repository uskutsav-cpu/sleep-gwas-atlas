#!/usr/bin/env python3
"""Rerun nine retained-input LDSC fits in a new directory; preserve saved fits."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import numpy as np

p = Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument('--output-dir', required=True)
args = ap.parse_args()
out = Path(args.output_dir).resolve()
out.mkdir(exist_ok=False)
root = Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas')
python = root/'.ldsc-env/bin/python3.9'
reference = Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr')
data = Path('/Volumes/Extreme SSD/brain6-physical-factor-20261101')
def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

r = json.loads((p/'prepare_receipt.json').read_text())
q = json.loads((p/'questionnaire_receipt.json').read_text())
expected = r['output_sha256'] | {'questionnaire': q['output_sha256']}
inputs = {str(data/(k+'.common.sumstats.gz')): v for k,v in expected.items()}
for f, v in inputs.items():
    assert sha(Path(f)) == v, ('input hash mismatch', f)
reference_files = [reference/(str(c)+suffix) for c in range(1,23)
                   for suffix in ('.l2.ldscore.gz', '.l2.M_5_50')]
code_files = sorted((root/'ldsc').rglob('*.py'))
assert code_files and len(reference_files) == 44
dependencies = {str(f): sha(f) for f in reference_files+code_files}
receipt = {'status': 'running', 'output_directory': str(out),
           'scope': 'Rerun LDSC from retained munged common-support inputs; not original raw GWAS ingestion or independent validation.',
           'input_sha256': inputs, 'reference_and_code_sha256': dependencies,
           'python_executable': str(python), 'fits': []}
receipt_path = out/'receipt.json'
def save():
    receipt_path.write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n')
save()
env = os.environ.copy()
env.update(MPLCONFIGDIR=str(out), OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
pairs = [('walking','grip'), ('walking','activity'), ('grip','activity')]
pairs += [(t,s) for s in ('insomnia','questionnaire') for t in ('walking','grip','activity')]
try:
    for a,b in pairs:
        name = a+'_'+b
        cmd = [str(python), str(root/'ldsc/ldsc.py'), '--rg',
               str(data/(a+'.common.sumstats.gz'))+','+str(data/(b+'.common.sumstats.gz')),
               '--ref-ld-chr', str(reference)+'/', '--w-ld-chr', str(reference)+'/',
               '--print-delete-vals', '--out', str(out/name)]
        with (out/(name+'.stdout.txt')).open('w') as log:
            subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        text = (out/(name+'.log')).read_text()
        original = (p/(name+'.log')).read_text()
        pattern = r'Genetic Correlation: ([\d.-]+) \(([\d.]+)\)'
        observed = re.search(pattern,text).groups()
        assert observed == re.search(pattern,original).groups(), name
        assert 'Analysis finished at' in text
        assert int(re.search(r'(\d+) SNPs with valid alleles',text).group(1)) == 1053614
        maximum = 0.
        for key in ('hsq1','hsq2','gencov'):
            suffix = a+'.common.sumstats.gz_'+b+'.common.sumstats.gz.'+key+'.delete'
            actual = np.loadtxt(out/(name+suffix))
            expected_values = np.loadtxt(p/(name+suffix))
            assert actual.shape == expected_values.shape == (200,)
            error = float(np.max(np.abs(actual-expected_values)))
            maximum = max(maximum,error)
            assert error < 1e-12, (name,key,error)
        receipt['fits'].append({'pair':name,'rg':float(observed[0]),'se':float(observed[1]),
                                'maximum_saved_delete_difference':maximum,'command':cmd})
        save()
        print(name, 'MATCH', flush=True)
    for f,v in inputs.items() | dependencies.items():
        assert sha(Path(f)) == v, ('dependency changed during replay',f)
    receipt['status'] = 'passed'
    receipt['all_nine_point_estimates_ses_and_27_delete_vectors_match'] = True
    save()
    print('PASSED', receipt_path, flush=True)
except BaseException as exc:
    receipt['status'] = 'failed'
    receipt['error'] = repr(exc)
    save()
    raise

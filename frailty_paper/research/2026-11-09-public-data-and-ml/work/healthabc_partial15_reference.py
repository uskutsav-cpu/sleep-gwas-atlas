#!/usr/bin/env python3
"""Compare Health ABC array subset to complete 15-marker weights in EUR reference genotypes."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

root = Path(__file__).resolve().parent.parent
score_path = root / 'outputs/unvalidated_15_variant_score_prototype.csv'
coverage_path = root / 'outputs/healthabc_geo_platform_coordinate_screen.json'
genotype_path = root / 'work/stage6_hfrs_sparse_score_EUR503_genotypes.npz'
score = list(csv.DictReader(score_path.open()))
screen = json.loads(coverage_path.read_text())
assert screen['complete_platform_table']
matched = {x['SNP_ID'] for x in screen['matches']}
with np.load(genotype_path, allow_pickle=True) as data:
    genotypes = np.asarray(data['genotypes'], dtype=float)
    ids = data['variants'].astype(str).tolist()
assert ids == [r['variant'] for r in score] and genotypes.shape == (503, 15) and np.isfinite(genotypes).all()
mask = np.array([r['variant'] in matched for r in score])
assert sum(mask) == 5
results = {}
for name in ('FI_beta_weight', 'FFS_beta_weight_inferred', 'insomnia_beta_weight'):
    weights = np.array([float(r[name]) for r in score])
    full = genotypes @ weights
    partial = genotypes[:, mask] @ weights[mask]
    results[name] = {'correlation': float(np.corrcoef(full, partial)[0, 1]),
                     'full_sd': float(full.std(ddof=1)), 'omitted_sd': float((full-partial).std(ddof=1))}
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
report = {'reference':'1000 Genomes Phase 3 EUR, 503 samples',
          'subset_rsids': sorted(r['variant'] for r in score if r['variant'] in matched),
          'results':results, 'input_sha256':{'score':digest(score_path),'coverage':digest(coverage_path),'genotypes':digest(genotype_path)},
          'limitations':'Reference genotype arithmetic only; no Health ABC data, target genotype QC, imputation, outcome, prediction or clinical validation. FFS effect orientation is inferred.'}
(root / 'outputs/healthabc_partial15_reference.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))

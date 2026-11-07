#!/usr/bin/env python3
"""Synthetic checks of summary subtraction; contains no participant data."""
import json
import math
from pathlib import Path
from statistics import NormalDist
import numpy as np

rng = np.random.default_rng(20261106)
n = 4000
g = rng.binomial(2, .3, n)
c = rng.normal(size=n)
shared = rng.normal(size=n)
items = np.column_stack([
    (.15*g + .2*c + .6*shared + rng.normal(size=n) > t).astype(float)
    for t in (-.6, -.2, .1, .3, .5)
])
full = items.sum(axis=1)
exhaustion = items[:, 0]
physical = full - exhaustion
design = np.column_stack([np.ones(n), g, c])
inverse = np.linalg.inv(design.T @ design)
ys = np.column_stack([full, exhaustion, physical])
betas = inverse @ design.T @ ys
residual = ys - design @ betas
sampling = inverse[1, 1] * (residual.T @ residual) / (n - design.shape[1])
reconstructed_variance = sampling[0, 0] + sampling[1, 1] - 2*sampling[0, 1]
assert np.max(abs(betas[:, 2] - (betas[:, 0]-betas[:, 1]))) < 1e-12
assert abs(sampling[2, 2]-reconstructed_variance) < 1e-12

def int_score(y):
    _, inverse_rank, counts = np.unique(y, return_inverse=True, return_counts=True)
    midrank = np.cumsum(counts) - (counts-1)/2
    z = np.array([NormalDist().inv_cdf((r-.5)/len(y)) for r in midrank])
    return z[inverse_rank]

int_betas = (inverse @ design.T @ np.column_stack([
    int_score(full), int_score(exhaustion), int_score(physical)]))[1]
assert abs(int_betas[2]-(int_betas[0]-int_betas[1])) > 1e-4

# Exchangeable physical indicators: conditional on count K, distribute ones
# uniformly across the four positions. Entries are (probability, K, E).
base = [(0.5, 0, 0), (0.25, 2, 1), (0.25, 4, 0)]
alternative = [(0.5, 0, 0), (0.25, 3, 1), (0.25, 3, 0)]
def moments(dist):
    return {
        'each_physical_item_probability': sum(p*k/4 for p,k,e in dist),
        'exhaustion_probability': sum(p*e for p,k,e in dist),
        'mean_five_item_count': sum(p*(k+e) for p,k,e in dist),
        'full_frailty_case_probability': sum(p for p,k,e in dist if k+e >= 3),
        'four_item_frailty_case_probability': sum(p for p,k,e in dist if k >= 3),
    }
a,b = moments(base),moments(alternative)
for key in a:
    if key != 'four_item_frailty_case_probability':
        assert a[key] == b[key]
assert a['four_item_frailty_case_probability'] == .25
assert b['four_item_frailty_case_probability'] == .5
result = {
    'synthetic_only': True, 'seed': 20261106, 'n': n,
    'linear_snp_beta_full_exhaustion_physical': betas[1].tolist(),
    'maximum_linear_coefficient_identity_error': float(np.max(abs(betas[:,2]-(betas[:,0]-betas[:,1])))),
    'direct_physical_beta_se': float(math.sqrt(sampling[2,2])),
    'covariance_aware_reconstructed_se': float(math.sqrt(reconstructed_variance)),
    'incorrect_independence_se': float(math.sqrt(sampling[0,0]+sampling[1,1])),
    'separately_int_snp_betas_full_exhaustion_physical': int_betas.tolist(),
    'incorrect_int_subtraction': float(int_betas[0]-int_betas[1]),
    'binary_counterexample_g0': a, 'binary_counterexample_g1': b,
    'binary_deleted_score_log_odds_ratio': math.log(3),
    'qualification': 'Mathematical feasibility check, not observed frailty genetics or a validation cohort.'
}
out = Path(__file__).with_name('linear_route_check.json')
out.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
print(json.dumps(result, indent=2))

import importlib.util
import math
from pathlib import Path
import numpy as np
from scipy.stats import norm

HERE=Path(__file__).resolve().parents[1]
def module(name):
    spec=importlib.util.spec_from_file_location(name,HERE/'scripts'/f'{name}.py')
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def test_covariance_bound_all_admissible_correlations():
    m=module('global_contrasts')
    for s1,s2,delta in [(.03,.025,.2),(.1,.1,0),(.02,.08,-.3)]:
        semax,pmax=m.bounded_contrast(delta,s1,s2)
        for c in np.linspace(-1,1,101):
            se=math.sqrt(max(0,s1*s1+s2*s2-2*c*s1*s2))
            p=2*norm.sf(abs(delta)/se) if se else (0 if delta else 1)
            assert se<=semax+1e-14 and p<=pmax+1e-14

def test_qc_rejects_low_h2_and_invalid_rg():
    m=module('global_contrasts')
    row=dict(rg=.4,se=.03,sleep_h2=.1,sleep_h2_se=.01,h2_disorder=.2,h2_disorder_se=.02)
    assert m.marginal_eligible(row)
    assert not m.marginal_eligible(dict(row,rg=1.2))
    assert not m.marginal_eligible(dict(row,sleep_h2=.001))
    assert not m.marginal_eligible(dict(row,se=0))

def test_coloc_distinct_variants_against_explicit_pair_products():
    m=module('molecular_sensitivity');l1=np.array([2.,.1,-1]);l2=np.array([-.4,1.,.7])
    a=np.exp(l1);b=np.exp(l2)
    weights=np.array([1,1e-4*sum(a),1e-4*sum(b),1e-8*sum(a[i]*b[j] for i in range(3) for j in range(3) if i!=j),1e-5*sum(a*b)])
    np.testing.assert_allclose(m.posterior(l1,l2,1e-4,1e-4,1e-5),weights/weights.sum(),atol=1e-14)

def test_coloc_sdy_known_regression_identity():
    m=module('molecular_sensitivity');maf=np.array([.1,.2,.3]);n=np.array([100,150,200]);sd=2.5
    v=sd**2/(2*n*maf*(1-maf))
    assert abs(m.estimate_sdy(v,maf,n)-sd)<1e-14

def test_logabf_sign_invariance_and_finite_large_evidence():
    m=module('molecular_sensitivity');beta=np.array([.1,-.3,.4]);se=np.array([.02,.08,.1])
    np.testing.assert_allclose(m.logabf(beta,se,.2),m.logabf(-beta,se,.2),atol=1e-14)
    p=m.posterior(np.array([700.,2.,0]),np.array([701.,1.,0]),1e-4,1e-4,1e-5)
    assert np.all(np.isfinite(p)) and abs(p.sum()-1)<1e-12

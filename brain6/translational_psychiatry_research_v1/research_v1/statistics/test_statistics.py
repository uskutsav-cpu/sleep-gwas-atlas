"""Regression tests for source identity, algebra, and inference counterexamples."""
import json, math
import numpy as np
import pandas as pd
from scipy.stats import norm
from calibration import HERE, native_fixture, reconstructed, exact_ld, verify_freeze
from tail_diagnostics import product_cumulants

def test_frozen_calibration_code_and_inputs_unchanged():
    assert verify_freeze()['no_simulation_outcomes_read']

def test_patched_source_matches_exact_historical_provenance():
    x=pd.read_csv(HERE/'historical_source_identity.tsv',sep='\t')
    assert len(x)==6
    assert set(x.status)=={'EXACT_MATCH'}

def test_independent_reconstruction_matches_exact_upstream_numerics():
    for ar in [.2,.6,.9,.97]:
        _,_,e,_=native_fixture(ar,.001,.2)
        assert e['rho']<1e-12 and e['var']<1e-15 and e['p']<1e-10

def test_scalar_n_rescaling_does_not_certify_model_semantics():
    baseline,_,_,(tx,ty,d,meanld,hc1,hc2)=native_fixture()
    e=reconstructed(tx[:1],ty[:1],d,386533*2,225534*2,np.array([hc1/2]),np.array([hc2/2]),meanld,.2)
    assert math.isclose(e['rho'][0],baseline.rho.iloc[0]/2,rel_tol=1e-12)
    assert math.isclose(e['var'][0],baseline['var'].iloc[0]/4,rel_tol=1e-12)
    assert abs(e['p'][0]-baseline.p.iloc[0])<1e-10

def test_moment_and_gls_variances_do_not_match_in_strong_ld():
    m=160; n1=386533; n2=225534; h=.01
    _,d,_,l2=exact_ld(m,.9)
    av=d+n1*h*d*d/m; bv=d+n2*h*d*d/m
    moment=np.sum(av*bv)/(n1*n2*np.mean(l2)**2)
    gls=m*m/(n1*n2*np.sum(d**4/(av*bv)))
    assert moment/gls>10
    # Weighted variance cannot be assigned to the unchanged unweighted moment point estimate.
    assert not math.isclose(moment,gls,rel_tol=.1)

def test_product_cumulant_formula_known_cases():
    v,k4=product_cumulants(1.,1.,0.); assert v==1. and k4==6.
    v,k4=product_cumulants(1.,1.,1.); assert v==2. and k4==48.

def test_unknown_sampling_covariance_contrast_bound_is_conservative():
    s1,s2=.03,.09; delta=.27
    maximum=2*norm.sf(abs(delta)/(s1+s2))
    for r in np.linspace(-1.,1.,101):
        var=s1*s1+s2*s2-2*r*s1*s2
        assert var<=(s1+s2)**2+1e-15
        assert 2*norm.sf(abs(delta)/math.sqrt(var))<=maximum+1e-15

def test_frozen_family_not_reduced_by_missing_blocks():
    x=json.loads((HERE/'empirical_diagnostics_summary.json').read_text())
    assert x['family_denominator']==8465
    assert math.isclose(x['alpha'],.05/8465)

def test_negative_determinant_does_not_create_bounded_correlation():
    x=pd.read_csv(HERE/'empirical_block_diagnostics.tsv',sep='\t')
    valid=x.both_h2_positive
    assert int((valid&(x.moment_parameter_matrix_determinant<0)).sum())==723
    assert int(x.derived_correlation_out_of_bounds.sum())==723
    assert not ((x['corr']==1)&(x.moment_parameter_matrix_determinant<0)).any()

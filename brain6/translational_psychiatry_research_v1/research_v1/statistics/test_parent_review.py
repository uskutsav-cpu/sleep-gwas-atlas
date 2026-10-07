"""Scientific-inference regression checks for the independently reviewed parent analyses."""
import json, math
import numpy as np
import pandas as pd
import pytest
from calibration import HERE, ROOT, sha
from review_parent_analyses import DATA, exclusive_pair_logsum, separate_logabf

def test_all_global_pairs_reviewed_with_original72_denominator():
    r=pd.read_csv(HERE/'parent_global_contrast_checks.tsv',sep='\t')
    assert len(r)==72 and r.reviewed_admitted.all() and r.reviewed_rejected.sum()==16
    assert r[[c for c in r if c.endswith('abs_error')]].to_numpy().max()<1e-12
    grid=pd.read_csv(DATA/'global_contrasts/covariance_sensitivity.tsv',sep='\t')
    observed=pd.read_csv(DATA/'global_contrasts/contrasts.tsv',sep='\t')
    joined=grid.merge(observed[['sleep_trait','disorder1','disorder2','maximum_normal_model_P']],on=['sleep_trait','disorder1','disorder2'])
    assert (joined.normal_model_P<=joined.maximum_normal_model_P+1e-12).all()

def test_explicit_distinct_variant_configuration_sum():
    a=np.array([2.,3.,5.]); b=np.array([7.,11.,13.])
    assert math.isclose(math.exp(exclusive_pair_logsum(np.log(a),np.log(b))),198.,rel_tol=1e-12)

def test_ABF_sign_invariance_cannot_validate_strand_harmonization():
    beta=np.array([.01,-.05,.2]); se=np.array([.01,.02,.03])
    assert np.array_equal(separate_logabf(beta,se,.2),separate_logabf(-beta,se,.2))

def test_original_molecular60_has_cross_language_numeric_agreement():
    py=pd.read_csv(HERE/'parent_molecular_configuration_checks.tsv',sep='\t')
    ref=pd.read_csv(HERE/'reference_R_posterior_checks.tsv',sep='\t')
    assert len(py)==len(ref)==60
    assert py.max_pp_abs_error.max()<1e-10
    assert ref.max_posterior_absolute_error.max()<1e-10
    assert not py.H4_ge08.any()

def test_palindromic_stress_preserves_all_contexts_and_input_identity():
    directory=HERE/'molecular_palindromic_diagnostic';cfg=json.loads((directory/'protocol.json').read_text())
    amendment=json.loads((directory/'storage_amendment.json').read_text())
    assert sha(directory/'frozen_analysis_code_v1.py.txt')==cfg['analysis_code_sha256']==amendment['original_code_sha256']
    assert sha(HERE/'molecular_palindromic_diagnostic.py')==amendment['corrected_code_sha256']
    for p,digest in cfg['dependency_code_sha256'].items():assert sha(HERE/p)==digest
    manifest=pd.read_csv(directory/'source_manifest.tsv',sep='\t');assert len(manifest)==4
    assert (manifest.retained_variants>=500).all()
    assert not list(directory.glob('*_nonpal.tsv.gz'))
    output=pd.read_csv(directory/'sensitivity_results.tsv',sep='\t')
    assert len(output)==60 and output.status.eq('DIAGNOSTIC_SINGLE_SIGNAL_ADMITTED').all()
    assert output.pp_h4.max()<.8

@pytest.mark.integration
def test_local_filtered_raw_slice_integration_requires_actual_inputs():
    directory=HERE/'molecular_palindromic_diagnostic'
    manifest=pd.read_csv(directory/'source_manifest.tsv',sep='\t');assert len(manifest)==4
    locations={r['file']:ROOT/r['local_storage_relative_path'] for r in json.loads((directory/'local_storage_manifest.json').read_text())['slices']}
    for r in manifest.itertuples():
        # Deliberately fail if actual required local data are missing; never silently skip.
        assert sha(locations[r.filtered_file])==r.filtered_sha256
        data=pd.read_csv(locations[r.filtered_file],sep='\t')
        pairs=data.variant.str.split('_').str[-2:].map(lambda a:frozenset(a))
        assert not pairs.isin([frozenset(['A','T']),frozenset(['C','G'])]).any()
        assert len(data)>=500 and data.gwas_p.min()<=5e-8 and data.qtl_p.min()<=5e-8

def test_independent_gaussian_quadratic_form_check_matches_parent_MC():
    table=pd.read_csv(HERE/'parent_independent_variance_checks.tsv',sep='\t')
    assert len(table)==2 and table.relative_error.max()<1e-12
    assert table.analytic_vs_parent_MC_in_MCSE.max()<3

def test_native_fine_mapping_derived_posteriors_and_credible_checks():
    summary=json.loads((HERE/'native_fine_mapping_review_summary.json').read_text())
    assert summary['native_fits_reviewed_from_RDS']==8
    assert summary['signal_pair_posterior_conditions_reviewed']==12
    assert summary['minimum_conditional_CS_coverage']>=.95
    assert summary['minimum_CS_purity']>=.5
    assert summary['minimum_ELBO_step']>=-1e-6
    assert summary['max_PIP_error']<1e-12 and summary['max_configuration_posterior_error']<1e-10
    assert summary['independent_two_trait_replications']==0
    posterior=pd.read_csv(HERE/'native_coloc_configuration_and_prior_checks.tsv',sep='\t')
    assert len(posterior)==12 and (posterior.all_hypotheses_sum-1).abs().max()<1e-12
    assert posterior.max_prior_reweighting_error.max()<1e-10

@pytest.mark.integration
def test_local_native_fit_and_LD_source_identity_integration():
    manifest=json.loads((HERE/'native_fine_mapping_review_input_manifest.json').read_text())
    assert len(manifest['local_inputs'])==13
    for row in manifest['local_inputs']:
        # Missing real fitted objects/BFs/LD causes failure when integration is requested.
        assert sha(ROOT/row['path_relative_to_repo'])==row['sha256']

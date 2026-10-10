"""Seal the narrow completed aggregate arithmetic audit; run no science builder."""
from pathlib import Path
import datetime
import hashlib
import json

P = Path(__file__).resolve().parents[1]
R = P/'reviews'
PREFIX = 'independent_shared_sleep_cov0_diagnostic41'
CHECKER = R/(PREFIX+'_checker_v1.py')
RECEIPT = R/(PREFIX+'_receipt_v1.json')
ARITHMETIC = R/(PREFIX+'_arithmetic_v1.tsv')
MD = R/(PREFIX+'_review_v1.md')
REPORT = R/(PREFIX+'_review_v1.json')
SEAL = R/(PREFIX+'_review_seal_v1.json')


def sha(path):
    assert path.is_file() and not path.is_symlink()
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path,value):
    with path.open('x') as stream:
        json.dump(value,stream,indent=2,allow_nan=False)
        stream.write('\n')


receipt=json.loads(RECEIPT.read_text())
assert receipt['all_checks_pass'] is True
assert all(sha(Path(q))==v for q,v in receipt['exact_current_input_builder_table_receipt_sha256'].items())
assert sha(ARITHMETIC)==receipt['output_arithmetic_table_sha256']
report=dict(schema='independent_shared_sleep_cov0_diagnostic41_review_v1',
    reviewed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    verdict='QUALIFIED_AGGREGATE_ARITHMETIC_AND_JOIN_PASS',
    material_numerical_join_or_qualification_errors=[],
    exact_current_input_builder_table_receipt_sha256=receipt['exact_current_input_builder_table_receipt_sha256'],
    independent_checker_sha256=sha(CHECKER),independent_arithmetic_receipt_sha256=sha(RECEIPT),
    independent_arithmetic_table_sha256=sha(ARITHMETIC),
    source_cardinalities=receipt['source_cardinalities'],source_rows_checked=41,
    all29_output_fields_checked=True,exact_join_membership_and_original_order_pass=True,
    native_marginal_inputs_not_recertified=True,
    numeric_method=receipt['numeric_method'],presentation_only_tolerance=receipt['arithmetic_check_tolerance'],
    max_absolute_error_by_field=receipt['max_absolute_error_by_field'],
    nominal_cov0_diagnostics_all41=16,nominal_cov0_diagnostics_qualified23=7,
    locked_qualified_class_count=23,locked_below_threshold_class_count=18,
    original7_annotation_preserved_and_equals_current_qualified_nominal7=True,
    extra9_nominal_diagnostics_in_locked_below_threshold18=True,
    shared_sleep=True,sampling_covariance_estimated=False,
    confirmed_heterogeneity=False,fully_independent_two_trait_replication=False,
    new_test_family_primary_claim_or_correction=False,arbitrary_covariance_grid=False,
    scientific_threshold_source_scope_changes=False,
    permitted_claim='These exact marginal-input values yield the reported41-row unverified Cov0 diagnostic arithmetic.',
    blocked_claims=['Confirmed heterogeneity or effect differences','Estimated or corrected sampling covariance',
        'Fully independent two-trait replication','New primary association or new biological hypothesis',
        'Causal explanation, mechanistic inference, novelty or power comparisons','New confirmatory41-test family'],
    variance_identity='Var(rg_external-rg_discovery)=SE_external^2+SE_discovery^2-2Cov(rg_external,rg_discovery);this table substitutesCov=0 without estimating it.',
    covariance_sign_or_true_magnitude_inferred=False,
    selection_qualification='Retains every existing41 estimate from the historical217 screen;23-class and7-membership are historical annotations. Displaying41 does not turn a selected historical screen into a new confirmatory family.',
    source_phenotype_ancestry_ascertainment_overlap_QC_limitations_not_resolved=True,
    native190_numerical_audits_or_old_suites_repeated=False,
    native_deletevectors_reference_raw_GWAS_runtime_body_reads=False,
    native_fits_workers_network_transfers_meters_censuses_mutex_calls=0,
    original_build_receipt_independent_review_pending_flag_preserved=True,
    new_report_is_additive_current_review_evidence=True,
    figures_reviewed_or_admitted=False,
    figure_requirements='Consume exact sealed41table/receipt identities, retain all41/23/7 annotations and explicit Cov0 unverified/sharedsleep labels; two figure pages need separate visual review.',
    remaining_human_methods_decision='A defensible sampling covariance estimator and its calibration/comparability must be independently validated before treating these differences as confirmed heterogeneity; no covariance-sign assumption or arbitrarygrid substitution admitted.',
    resource_observations=dict(elapsed_seconds=receipt['elapsed_seconds'],max_RSS_bytes=receipt['max_RSS_bytes']),
)
save(REPORT,report)
builder=P/'scripts/106_build_shared_sleep_zero_cov_diagnostic_v4_1.py'
table=P/'tables/shared_sleep_cov0_diagnostic41_v4_1.tsv'
build_receipt=P/'logs/shared_sleep_cov0_diagnostic41_build_receipt_v4_1.json'
master=P.parent/'discovery_extension/sleep_submission_evidence_v1/tables/HETEROGENEITY_MASTER.tsv'
seven=P.parent/'discovery_extension/sleep_submission_evidence_v1/tables/heterogeneity_7.tsv'
text=f"""# Independent all41 zero-covariance diagnostic review

**QUALIFIED_AGGREGATE_ARITHMETIC_AND_JOIN_PASS.** No numerical, join or qualification error was found in the exact new [builder]({builder}), [41-row table]({table}) and [build receipt]({build_receipt}). This checks aggregate-table arithmetic and annotation only. It does not repeat or re-certify native190 audits, source compatibility, sampling covariance, estimator calibration or clinical semantics. No fits, native vectors, GWAS/reference/runtime bodies, locks, transfers or SSD/resource census were used.

The table SHA is `0698b9493637176c03b975bc36e9fbd6e2c5baf6dfdc545486c004df693149bb`, build receipt `c2f623a1f4759a68c8de4842b7a29f8ca524873048aaa5bab2e5c63c68c21ca4`, builder `72a8bb03b2cb10a774d3365a3d1c8f3e25cf896ecdb84868d79447bb3e823a1a`. The [independent checker]({CHECKER}), [receipt]({RECEIPT}) and [per-row decimal/line evidence]({ARITHMETIC}) bind all four exact aggregate sources and all41 rows, without importing or rerunning builder106.

The discovery aggregate has1200 distinct sleep/outcome keys. Native validation, [historical master]({master}) and new diagnostic table each contain41 unique pair IDs with identical membership. Each diagnostic row joins the matching sleep trait, original extension outcome and external source; historical phenotype/classification, native ratio rg and marginal SE, and source row identity all agree. The original23 qualified external outcome-side validations and18 directionally concordant below-threshold pairs remain unchanged. Output ordering matches the declared sleep/phenotype/pair sort. All29 columns were checked, including15 numerical fields, seven Boolean qualifications and exact identity/interpretation strings. The seven full [historical rows]({seven}) agree with their master records and their membership is retained; historical covariance grids/correction columns were not used by the new builder.

For each row, independently compute D=external rg-discovery rg, SE0=sqrt(SEexternal²+SEdiscovery²), Z0=D/SE0 and CI0=D±1.959963984540054×SE0. An80-digit Decimal calculation starting from literal aggregate input strings verifies these, both historical-printing differences, and P0 using an independent erf power series (1-erf(|Z0|/sqrt2)), rather than builder math.hypot/math.erfc. All328 derived values pass the declared presentation-only arithmetic tolerance. Maximum absolute discrepancy is8.76×10^-16 for Z and2.60×10^-16 for P; this is ordinary floating-point representation/operation scale, not a scientific reproduction threshold. Maximum |Z| is4.488, within the series's bounded numerical domain. The closest displayed P to0.05 remains0.001072 away, so the nominal flags are unaffected by these errors.

There are16 unverified nominal Cov0 diagnostics among all41, consisting of seven among the historical qualified23 and nine among the below-threshold18. The seven qualified nominal members are exactly the original seven-row annotation. These counts describe the frozen selected estimates; they do not confirm effect differences, show a new discovery, constitute a corrected41-test family or establish biological contrasts. All41 are retained without a new filter or correction. The build receipt's review-pending flag is preserved as build-time history; this additive seal records the current independent arithmetic review.

The actual variance identity is Var(external rg-discovery rg)=SEexternal²+SEdiscovery²-2Cov(external rg,discovery rg). The displayed SE, CI and P substitute0 for this unestimated cross-estimator covariance. Reusing the sleep GWAS leaves this assumption unverified and prevents a fully independent two-trait replication interpretation. It does not establish the sign or magnitude of the true covariance, so no claim that Cov0 is conservative or anti-conservative follows. Exact arithmetic cannot resolve source/phenotype/ancestry/ascertainment comparability, participant overlap, shared-source uncertainty, historical source/QC limitations or selected-screen/winner's-curse concerns.

All rows correctly label `shared_sleep_source=True`, `sampling_covariance_established=False`, `confirmed_effect_difference=False`, `fully_independent_two_trait_replication=False`, `new_test_family_or_primary_claim=False` and `UNVERIFIED_ZERO_COVARIANCE_DIAGNOSTIC_ONLY`; every derived interval/P column names its unvalidated Cov0 condition. The original1200-test and217-candidate families and thresholds are not modified. No arbitrary covariance grid is promoted to evidence. A justified and independently calibrated sampling-covariance estimator, together with compatible sources/estimands and an appropriate inference protocol, is still required before any confirmed heterogeneity claim.

The proposed two diagnostic figure pages are not reviewed by this seal. They should consume the exact sealed table and receipt hashes, retain the41/23/7 distinctions, state that shared covariance is unestimated, and receive separate visual verification. No causal, mechanistic, novelty, power or new biological conclusion is admitted by this table review.
"""
with MD.open('x') as stream:stream.write(text)
files=dict(receipt['exact_current_input_builder_table_receipt_sha256'])
for q in [CHECKER,RECEIPT,ARITHMETIC,REPORT,MD,Path(__file__)]:files[str(q)]=sha(q)
assert all(sha(Path(q))==v for q,v in files.items())
save(SEAL,dict(schema='independent_shared_sleep_cov0_diagnostic41_review_seal_v1',
    verdict=report['verdict'],file_sha256=files,regular_file_count=len(files),
    native190_or_old_suite_reruns=False,raw_vector_reference_runtime_body_reads=False,
    covariance_validated=False,heterogeneity_confirmed=False,figure_review=False))
assert all(sha(Path(q))==v for q,v in files.items())
print(json.dumps(dict(verdict=report['verdict'],md=str(MD),md_sha256=sha(MD),
    report=str(REPORT),report_sha256=sha(REPORT),seal=str(SEAL),seal_sha256=sha(SEAL),
    file_count=len(files))))

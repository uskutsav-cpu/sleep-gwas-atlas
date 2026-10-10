"""Independent all41 aggregate-table join and decimal arithmetic check.

Does not import/run builder106, read native captures/deletevectors/GWAS/reference
bodies, perform native audits/fits, or acquire locks. Decimal erf uses its power
series, independently of builder math.erfc/hypot. All numerical tolerances here
are presentation-arithmetic checks, not scientific admission thresholds.
"""
import ast
import csv
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
import resource
import time

P = Path(__file__).resolve().parents[1]
ROOT = P.parent
R = P / 'reviews'
started = time.monotonic()
FILES = {
    P/'tables/extension_native_full_precision_rg_v4.tsv': 'c7551a1d3e246234289865469729618e61bdb6c35c19cffd8b5946a61a06b85b',
    P/'tables/validation_native_full_precision_rg_v4.tsv': '794dc2f951268d673d667e8e5d2cab1558fc9cf34e06673e5e2fe8e2979ee65b',
    ROOT/'discovery_extension/sleep_submission_evidence_v1/tables/HETEROGENEITY_MASTER.tsv': '062cee7515f710da1003e6ab29b1b51b64bb23ce86e3c86b8135630073777810',
    ROOT/'discovery_extension/sleep_submission_evidence_v1/tables/heterogeneity_7.tsv': '900243199959a657a197e65a555ae62cec851504ccb5f512e29b80255867a363',
    P/'tables/shared_sleep_cov0_diagnostic41_v4_1.tsv': '0698b9493637176c03b975bc36e9fbd6e2c5baf6dfdc545486c004df693149bb',
    P/'logs/shared_sleep_cov0_diagnostic41_build_receipt_v4_1.json': 'c2f623a1f4759a68c8de4842b7a29f8ca524873048aaa5bab2e5c63c68c21ca4',
    P/'scripts/106_build_shared_sleep_zero_cov_diagnostic_v4_1.py': '72a8bb03b2cb10a774d3365a3d1c8f3e25cf896ecdb84868d79447bb3e823a1a',
}
PI = Decimal('3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679')
Z95 = Decimal('1.959963984540054')
OUT = R/'independent_shared_sleep_cov0_diagnostic41_receipt_v1.json'
ROWS_OUT = R/'independent_shared_sleep_cov0_diagnostic41_arithmetic_v1.tsv'


def bound():
    assert time.monotonic()-started < 60
    assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss < (64 << 20)


def sha(path):
    assert path.is_file() and not path.is_symlink()
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        rows = list(reader)
        assert all(None not in r and all(v is not None for v in r.values()) for r in rows)
        return reader.fieldnames, rows


def erfc_decimal(x):
    """At observed x<=3.174,80-digit alternating erf series is well conditioned."""
    assert Decimal(0) <= x < Decimal('3.2')
    power_term = x
    total = x
    for n in range(1, 500):
        power_term *= -(x*x)/Decimal(n)
        term = power_term/Decimal(2*n+1)
        total += term
        if n > 25 and abs(term) < Decimal('1e-72'):
            value = 1-2*total/PI.sqrt()
            assert 0 < value <= 1
            return value, n
    raise AssertionError('DECIMAL_SERIES_DID_NOT_CONVERGE')


def compare(observed, expected, is_p=False):
    value = Decimal(observed)
    assert value.is_finite() and expected.is_finite()
    # Decimal literals from aggregate TSVs are the independent reference.
    # Double rounding/operation errors need allowance; these do not change QC.
    tolerance = Decimal('5e-14')*abs(expected) + (Decimal('1e-24') if is_p else Decimal('2e-15'))
    error = abs(value-expected)
    assert error <= tolerance, (observed, str(expected), str(error), str(tolerance))
    return error


assert all(sha(q) == expected for q, expected in FILES.items())
_, ext = read(P/'tables/extension_native_full_precision_rg_v4.tsv')
_, val = read(P/'tables/validation_native_full_precision_rg_v4.tsv')
_, historical = read(ROOT/'discovery_extension/sleep_submission_evidence_v1/tables/HETEROGENEITY_MASTER.tsv')
_, seven = read(ROOT/'discovery_extension/sleep_submission_evidence_v1/tables/heterogeneity_7.tsv')
columns, diagnostic = read(P/'tables/shared_sleep_cov0_diagnostic41_v4_1.tsv')
assert len(ext) == 1200 and len(val) == len(historical) == len(diagnostic) == 41 and len(seven) == 7
extmap = {(x['sleep_trait'],x['outcome_trait']): (i+2,x) for i,x in enumerate(ext)}
valmap = {x['original_pair_id']: (i+2,x) for i,x in enumerate(val)}
oldmap = {x['pair_id']: (i+2,x) for i,x in enumerate(historical)}
sevenmap = {x['pair_id']: (i+2,x) for i,x in enumerate(seven)}
outmap = {x['pair_id']: (i+2,x) for i,x in enumerate(diagnostic)}
assert len(extmap) == 1200 and len(valmap) == len(oldmap) == len(outmap) == 41 and len(sevenmap) == 7
assert set(outmap) == set(valmap) == set(oldmap) and set(sevenmap) <= set(oldmap)
assert all(x == oldmap[k][1] for k,(_,x) in sevenmap.items())
expected_order = sorted(val, key=lambda x:(x['sleep_trait'],oldmap[x['original_pair_id']][1]['external_phenotype_name'],x['original_pair_id']))
assert [r['pair_id'] for r in diagnostic] == [r['original_pair_id'] for r in expected_order]
assert set(columns) == {
    'pair_id','sleep_trait','discovery_outcome_id','external_source_id','phenotype','locked_validation_class',
    'discovery_native_rg','discovery_native_SE','external_native_rg','external_native_SE',
    'difference_external_minus_discovery','assumed_sampling_covariance_for_diagnostic',
    'difference_SE_cov0_unvalidated','difference_CI95_low_cov0_unvalidated','difference_CI95_high_cov0_unvalidated',
    'difference_Z_cov0_unvalidated','difference_P_cov0_unvalidated','historical_printed_difference',
    'historical_printed_difference_SE_cov0','native_minus_historical_printed_difference',
    'native_minus_historical_printed_difference_SE_cov0','original_nominal7_membership',
    'current_cov0_nominal_P_lt_0_05_for_diagnostic_only','shared_sleep_source','sampling_covariance_established',
    'confirmed_effect_difference','fully_independent_two_trait_replication','new_test_family_or_primary_claim','interpretation',
}
records = []
max_errors = {}
with localcontext() as ctx:
    ctx.prec = 80
    for rowline,t in enumerate(diagnostic, 2):
        bound()
        pair = t['pair_id']
        vl,v = valmap[pair]
        hl,h = oldmap[pair]
        dl,d = extmap[(v['sleep_trait'],v['original_extension_trait_id'])]
        assert pair == v['sleep_trait']+'__'+v['original_extension_trait_id']
        assert h['sleep_trait'] == v['sleep_trait'] == v['original_sleep_trait'] == t['sleep_trait']
        assert h['extension_trait_id'] == v['original_extension_trait_id'] == d['outcome_trait'] == t['discovery_outcome_id']
        assert h['replication_source_id'] == v['outcome_trait'] == v['original_replication_source_id'] == t['external_source_id']
        assert t['phenotype'] == h['external_phenotype_name']
        assert t['locked_validation_class'] == v['historical_classification']
        assert v['sampling_covariance_corrected'] == v['independent_two_trait_replication'] == 'False'
        assert h['sleep_GWAS_reused'] == 'True' and h['heterogeneity_assumption'] == 'ZERO_COVARIANCE_ASSUMPTION_UNVERIFIED'
        assert t['original_nominal7_membership'] == str(pair in sevenmap)
        assert t['shared_sleep_source'] == 'True'
        for k in ['sampling_covariance_established','confirmed_effect_difference','fully_independent_two_trait_replication','new_test_family_or_primary_claim']:
            assert t[k] == 'False'
        assert t['interpretation'] == 'UNVERIFIED_ZERO_COVARIANCE_DIAGNOSTIC_ONLY'
        assert t['assumed_sampling_covariance_for_diagnostic'] == '0'
        for key,source,field in [('discovery_native_rg',d,'rg'),('discovery_native_SE',d,'se'),('external_native_rg',v,'rg'),('external_native_SE',v,'se')]:
            assert t[key] == source[field]
        for key,oldfield in [('historical_printed_difference','difference_replication_minus_discovery'),('historical_printed_difference_SE_cov0','difference_se_cov0')]:
            assert t[key] == h[oldfield]
        dr,ds,vr,vs = [Decimal(x) for x in (d['rg'],d['se'],v['rg'],v['se'])]
        assert all(x.is_finite() for x in (dr,ds,vr,vs)) and ds > 0 and vs > 0
        difference = vr-dr
        se = (vs*vs+ds*ds).sqrt()
        z = difference/se
        p,iterations = erfc_decimal(abs(z)/Decimal(2).sqrt())
        arithmetic = {
            'difference_external_minus_discovery': difference,
            'difference_SE_cov0_unvalidated': se,
            'difference_CI95_low_cov0_unvalidated': difference-Z95*se,
            'difference_CI95_high_cov0_unvalidated': difference+Z95*se,
            'difference_Z_cov0_unvalidated': z,
            'difference_P_cov0_unvalidated': p,
            'native_minus_historical_printed_difference': difference-Decimal(h['difference_replication_minus_discovery']),
            'native_minus_historical_printed_difference_SE_cov0': se-Decimal(h['difference_se_cov0']),
        }
        errors = {k:compare(t[k],x,k=='difference_P_cov0_unvalidated') for k,x in arithmetic.items()}
        for k,x in errors.items():max_errors[k] = max(max_errors.get(k,Decimal(0)),x)
        nominal = p < Decimal('0.05')
        assert t['current_cov0_nominal_P_lt_0_05_for_diagnostic_only'] == str(nominal)
        records.append(dict(pair_id=pair,table_row=rowline,validation_row=vl,discovery_row=dl,
            historical_master_row=hl,historical7_row=sevenmap[pair][0] if pair in sevenmap else '',
            locked_validation_class=t['locked_validation_class'],decimal_D=str(difference),decimal_SE0=str(se),
            decimal_CI95_low=str(arithmetic['difference_CI95_low_cov0_unvalidated']),
            decimal_CI95_high=str(arithmetic['difference_CI95_high_cov0_unvalidated']),decimal_Z0=str(z),decimal_P0=str(p),
            decimal_P_series_iterations=iterations,original_nominal7=pair in sevenmap,
            diagnostic_nominal=nominal,max_numeric_absolute_error=str(max(errors.values())),
            every_numeric_join_boolean_and_qualification_pass=True))
assert sum(r['diagnostic_nominal'] for r in records) == 16
assert sum(r['diagnostic_nominal'] and r['locked_validation_class']=='QUALIFIED_EXTERNAL_OUTCOME_SIDE_VALIDATION' for r in records) == 7
assert sum(r['locked_validation_class']=='QUALIFIED_EXTERNAL_OUTCOME_SIDE_VALIDATION' for r in records) == 23
assert sum(r['locked_validation_class']=='DIRECTIONALLY_CONCORDANT_BELOW_FROZEN_THRESHOLD' for r in records) == 18
assert {r['pair_id'] for r in records if r['original_nominal7']} == {r['pair_id'] for r in records if r['diagnostic_nominal'] and r['locked_validation_class']=='QUALIFIED_EXTERNAL_OUTCOME_SIDE_VALIDATION'}
receipt = json.loads((P/'logs/shared_sleep_cov0_diagnostic41_build_receipt_v4_1.json').read_text())
assert receipt['input_sha256'] == {str(q):v for q,v in list(FILES.items())[:4]}
assert receipt['script_sha256'] == FILES[P/'scripts/106_build_shared_sleep_zero_cov_diagnostic_v4_1.py']
assert receipt['table_sha256'] == FILES[P/'tables/shared_sleep_cov0_diagnostic41_v4_1.tsv']
assert receipt['table_path'] == str(P/'tables/shared_sleep_cov0_diagnostic41_v4_1.tsv')
for key,n in [('exact_locked_estimated_pairs',41),('original_qualified_outcome_positive_pairs',23),('original_nominal7_membership_preserved',7),('current_cov0_nominal_P_lt_0_05_all41',16),('current_cov0_nominal_P_lt_0_05_qualified23',7),('estimator_calls',0)]:
    assert receipt[key] == n
for key in ['sampling_covariance_established','confirmed_heterogeneity','new_BH_or_other_correction','arbitrary_covariance_grid','original_thresholds_or_membership_changed','GWAS_body_reads']:
    assert receipt[key] is False
assert receipt['independent_numerical_review_pending'] is True  # Historical build-time flag preserved.
assert receipt['max_abs_native_minus_historical_printed_difference'] == max(abs(float(t['native_minus_historical_printed_difference'])) for t in diagnostic)
assert receipt['max_abs_native_minus_historical_printed_difference_SE_cov0'] == max(abs(float(t['native_minus_historical_printed_difference_SE_cov0'])) for t in diagnostic)
tree = ast.parse((P/'scripts/106_build_shared_sleep_zero_cov_diagnostic_v4_1.py').read_text())
strings = {n.value for n in ast.walk(tree) if isinstance(n,ast.Constant) and isinstance(n.value,str)}
assert not any('heterogeneity_p_rho_' in x or 'heterogeneity_bh_41' in x or 'heterogeneity_bonferroni_41' in x for x in strings)
assert all(sha(q) == v for q,v in FILES.items())
with ROWS_OUT.open('x',newline='') as stream:
    writer = csv.DictWriter(stream,fieldnames=list(records[0]),delimiter='\t')
    writer.writeheader();writer.writerows(records)
result = dict(schema='independent_all41_cov0_aggregate_decimal_review_v1',all_checks_pass=True,
    exact_current_input_builder_table_receipt_sha256={str(q):v for q,v in FILES.items()},
    source_cardinalities=dict(native_discovery=1200,native_validation=41,historical_master=41,historical7=7,diagnostic=41),
    unique_complete_join_order_historical_membership_every_field_qualification=True,
    numeric_method='80-digit Decimal from literal aggregate input strings; sqrt(sum of squares); independent erf power series/1-erf, not builder math.erfc/hypot;80digit context/1e-72 termination.',
    arithmetic_check_tolerance='abs_error<=5e-14*abs(expected)+2e-15 (P uses1e-24 absolute term);presentation only,no scientific threshold changes.',
    max_absolute_error_by_field={k:str(v) for k,v in max_errors.items()},
    numerical_fields_per_row=15,derived_decimal_fields_per_row=8,
    qualified_count=23,below_frozen_threshold_count=18,nominal_cov0_all41=16,nominal_cov0_qualified23=7,
    original7_equals_current_qualified_nominal7=True,
    nominal_below_frozen_threshold18=9,
    minimum_distance_to_nominal0p05=str(min(abs(Decimal(t['difference_P_cov0_unvalidated'])-Decimal('0.05')) for t in diagnostic)),
    assumed_covariance_not_estimated=True,confirmed_heterogeneity=False,new_test_family_or_correction=False,
    native190_audits_or_old_suites_rerun=False,GWAS_nativevector_reference_runtime_body_reads=False,
    nativefit_network_transfer_mutex_operations=0,
    builder_receipt_pending_flag_preserved_as_build_time_history=True,
    output_arithmetic_table=str(ROWS_OUT),output_arithmetic_table_sha256=sha(ROWS_OUT),
    elapsed_seconds=time.monotonic()-started,max_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    resource_bounds=dict(seconds=60,RSS_bytes=64 << 20,workers=0),
    limitations='Checks arithmetic, joins and exact qualification of marginal aggregate inputs only. Does not re-certify native190, source semantics, overlap, sampling covariance or estimator calibration. Shared sleep means zero cross-estimator covariance is unverified; no direction or magnitude of true covariance is inferred. Historical23 were outcome-side validations under the original217 family, not a newly selected confirmatory41-test family.')
with OUT.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
bound()
print(json.dumps(dict(all_checks_pass=True,rows=41,nominal_all41=16,nominal_qualified23=7,
    max_error_by_field=result['max_absolute_error_by_field'],elapsed_seconds=result['elapsed_seconds'],max_RSS_bytes=result['max_RSS_bytes'])))

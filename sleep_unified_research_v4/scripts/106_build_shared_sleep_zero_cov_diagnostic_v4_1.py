"""Describe all 41 locked effect differences under the unverified Cov=0 assumption.

This diagnostic cannot establish sampling covariance or confirmed heterogeneity.
No new test family, selection, estimator fit, or biological claim is introduced.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
from datetime import datetime, timezone

P = Path(__file__).resolve().parents[1]
ROOT = P.parent
INPUTS = {
    P/'tables/extension_native_full_precision_rg_v4.tsv': 'c7551a1d3e246234289865469729618e61bdb6c35c19cffd8b5946a61a06b85b',
    P/'tables/validation_native_full_precision_rg_v4.tsv': '794dc2f951268d673d667e8e5d2cab1558fc9cf34e06673e5e2fe8e2979ee65b',
    ROOT/'discovery_extension/sleep_submission_evidence_v1/tables/HETEROGENEITY_MASTER.tsv': '062cee7515f710da1003e6ab29b1b51b64bb23ce86e3c86b8135630073777810',
    ROOT/'discovery_extension/sleep_submission_evidence_v1/tables/heterogeneity_7.tsv': '900243199959a657a197e65a555ae62cec851504ccb5f512e29b80255867a363',
}
TABLE = P/'tables/shared_sleep_cov0_diagnostic41_v4_1.tsv'
RECEIPT = P/'logs/shared_sleep_cov0_diagnostic41_build_receipt_v4_1.json'


def sha(q):
    return hashlib.sha256(q.read_bytes()).hexdigest()


def main():
    assert not TABLE.exists() and not RECEIPT.exists()
    rows = {}
    for q, expected in INPUTS.items():
        assert q.is_file() and not q.is_symlink() and sha(q) == expected
        with q.open(newline='') as stream:
            rows[q.name] = list(csv.DictReader(stream, delimiter='\t'))
    discovery = {(r['sleep_trait'], r['outcome_trait']): r for r in rows['extension_native_full_precision_rg_v4.tsv']}
    validation = rows['validation_native_full_precision_rg_v4.tsv']
    historical = {r['pair_id']: r for r in rows['HETEROGENEITY_MASTER.tsv']}
    historical7 = {r['pair_id'] for r in rows['heterogeneity_7.tsv']}
    assert len(discovery) == 1200 and len(validation) == len(historical) == 41 and len(historical7) == 7
    assert {r['original_pair_id'] for r in validation} == set(historical)
    assert historical7 <= set(historical)
    result = []
    for v in sorted(validation, key=lambda x: (x['sleep_trait'], historical[x['original_pair_id']]['external_phenotype_name'], x['original_pair_id'])):
        pair = v['original_pair_id']
        old = historical[pair]
        d = discovery[(v['sleep_trait'], v['original_extension_trait_id'])]
        assert old['sleep_trait'] == v['sleep_trait'] == v['original_sleep_trait']
        assert old['extension_trait_id'] == v['original_extension_trait_id']
        assert old['replication_source_id'] == v['outcome_trait'] == v['original_replication_source_id']
        assert v['sampling_covariance_corrected'] == 'False' and v['independent_two_trait_replication'] == 'False'
        dr, ds, vr, vs = (float(d['rg']), float(d['se']), float(v['rg']), float(v['se']))
        assert all(math.isfinite(x) for x in (dr, ds, vr, vs)) and ds > 0 and vs > 0
        difference = vr - dr
        se0 = math.hypot(ds, vs)
        z0 = difference / se0
        p0 = math.erfc(abs(z0) / math.sqrt(2))
        old_difference = float(old['difference_replication_minus_discovery'])
        old_se0 = float(old['difference_se_cov0'])
        zcrit = 1.959963984540054
        result.append(dict(
            pair_id=pair, sleep_trait=v['sleep_trait'], discovery_outcome_id=v['original_extension_trait_id'],
            external_source_id=v['outcome_trait'], phenotype=old['external_phenotype_name'],
            locked_validation_class=v['historical_classification'],
            discovery_native_rg=dr, discovery_native_SE=ds,
            external_native_rg=vr, external_native_SE=vs,
            difference_external_minus_discovery=difference,
            assumed_sampling_covariance_for_diagnostic=0,
            difference_SE_cov0_unvalidated=se0,
            difference_CI95_low_cov0_unvalidated=difference-zcrit*se0,
            difference_CI95_high_cov0_unvalidated=difference+zcrit*se0,
            difference_Z_cov0_unvalidated=z0,
            difference_P_cov0_unvalidated=p0,
            historical_printed_difference=old_difference,
            historical_printed_difference_SE_cov0=old_se0,
            native_minus_historical_printed_difference=difference-old_difference,
            native_minus_historical_printed_difference_SE_cov0=se0-old_se0,
            original_nominal7_membership=pair in historical7,
            current_cov0_nominal_P_lt_0_05_for_diagnostic_only=p0<0.05,
            shared_sleep_source=True, sampling_covariance_established=False,
            confirmed_effect_difference=False, fully_independent_two_trait_replication=False,
            new_test_family_or_primary_claim=False,
            interpretation='UNVERIFIED_ZERO_COVARIANCE_DIAGNOSTIC_ONLY',
        ))
    assert len(result) == len({r['pair_id'] for r in result}) == 41
    assert sum(r['locked_validation_class'] == 'QUALIFIED_EXTERNAL_OUTCOME_SIDE_VALIDATION' for r in result) == 23
    for q, expected in INPUTS.items():
        assert sha(q) == expected
    with TABLE.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(result[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(result)
    receipt = dict(
        schema='all41_shared_sleep_cov0_descriptive_diagnostic_v1',
        recorded_utc=datetime.now(timezone.utc).isoformat(),
        input_sha256={str(q): digest for q, digest in INPUTS.items()},
        script_sha256=sha(Path(__file__)), table_path=str(TABLE), table_sha256=sha(TABLE),
        exact_locked_estimated_pairs=41, original_qualified_outcome_positive_pairs=23,
        original_nominal7_membership_preserved=7,
        current_cov0_nominal_P_lt_0_05_all41=sum(r['current_cov0_nominal_P_lt_0_05_for_diagnostic_only'] for r in result),
        current_cov0_nominal_P_lt_0_05_qualified23=sum(r['current_cov0_nominal_P_lt_0_05_for_diagnostic_only'] for r in result if r['locked_validation_class']=='QUALIFIED_EXTERNAL_OUTCOME_SIDE_VALIDATION'),
        max_abs_native_minus_historical_printed_difference=max(abs(r['native_minus_historical_printed_difference']) for r in result),
        max_abs_native_minus_historical_printed_difference_SE_cov0=max(abs(r['native_minus_historical_printed_difference_SE_cov0']) for r in result),
        formula='D=rg_external-rg_discovery; SE0=hypot(SE_external,SE_discovery); CI0=D +/- 1.959963984540054*SE0; P0=erfc(abs(D/SE0)/sqrt(2))',
        sampling_covariance_established=False, confirmed_heterogeneity=False,
        new_BH_or_other_correction=False, arbitrary_covariance_grid=False,
        original_thresholds_or_membership_changed=False, GWAS_body_reads=False, estimator_calls=0,
        scientific_qualification='All41 existing estimated pairs are retained. The sleep GWAS is reused and Cov0 is unverified, so these marginal-input effect-difference intervals and P values are descriptive diagnostics only. Historical7 is a preserved historical annotation. No result validates independence, corrected covariance, heterogeneity, causal mechanism or a new biological hypothesis.',
        independent_numerical_review_pending=True,
    )
    with RECEIPT.open('x') as stream:
        stream.write(json.dumps(receipt, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status':'DIAGNOSTIC41_BUILT_NOT_COVARIANCE_VALIDATED',
        'table_sha256':sha(TABLE), 'receipt_sha256':sha(RECEIPT),
        'current_cov0_nominal_all41':receipt['current_cov0_nominal_P_lt_0_05_all41'],
        'current_cov0_nominal_qualified23':receipt['current_cov0_nominal_P_lt_0_05_qualified23'],
        'confirmed_heterogeneity':False}))


if __name__ == '__main__':
    main()

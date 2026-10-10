#!/usr/bin/env python3
"""Read-only, portable checks of released aggregate tables; no SSD access."""
import collections
import csv
import hashlib
import json
import math
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1]
TABLES = {
    'core_native_full_precision_rg_v4.tsv': ('1928f02fdf9a56b37afdd67f7fb117dcae198098994236452d50cde04b355c77', 396),
    'extension_native_full_precision_rg_v4.tsv': ('c7551a1d3e246234289865469729618e61bdb6c35c19cffd8b5946a61a06b85b', 1200),
    'validation_native_full_precision_rg_v4.tsv': ('794dc2f951268d673d667e8e5d2cab1558fc9cf34e06673e5e2fe8e2979ee65b', 41),
    'independent_validation_admissibility217_v4_1.tsv': ('60592d55c3358587b1143fd1998c0ce5d89317b831e1f0da4e62bf8fa6f606dd', 217),
    'qualified_validation_prior_art23_v4_1.tsv': ('ea26c99bdda5361994911feb3ab84b8490e89ebefae5291a414144c84b740d8f', 23),
    'independent_sensitivity_62_full_precision_with_CI_v4_6.tsv': ('507c56848e98c9b87b7f7f83c4a924d44b5f8d1e764e60d12566afffc8695e61', 62),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    rows = {}
    observed = {}
    for name, (expected, count) in TABLES.items():
        path = PACKAGE / 'tables' / name
        require(path.is_file() and not path.is_symlink(), 'REGULAR_AGGREGATE_TABLE_REQUIRED')
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        require(digest == expected, 'EXACT_AGGREGATE_TABLE_IDENTITY_REQUIRED: ' + name)
        with path.open(newline='') as stream:
            rows[name] = list(csv.DictReader(stream, delimiter='\t'))
        require(len(rows[name]) == count, 'LOCKED_TABLE_CARDINALITY_REQUIRED: ' + name)
        observed[name] = digest
    for name in list(TABLES)[:3]:
        records = rows[name]
        require(len({(r['sleep_trait'], r['outcome_trait']) for r in records}) == len(records), 'UNIQUE_NATIVE_PAIR_KEYS_REQUIRED')
        for row in records:
            require(row['reproduction_status'] == 'PRINTED_PRECISION_CONCORDANT', 'NATIVE_PRINTED_PRECISION_STATUS_CHANGED')
            require(row['historical_full_precision_available'] == 'False', 'HISTORICAL_PRECISION_QUALIFICATION_CHANGED')
            rg, se, lower, upper = (float(row[key]) for key in ['rg', 'se', 'ci_lower_95', 'ci_upper_95'])
            require(all(math.isfinite(value) for value in [rg, se, lower, upper]) and se > 0 and lower < upper, 'FINITE_MARGINAL_ESTIMATE_AND_INTERVAL_REQUIRED')
    core = rows['core_native_full_precision_rg_v4.tsv']
    extension = rows['extension_native_full_precision_rg_v4.tsv']
    require(sum(r['native_fdr_pass_0_05'] == 'True' for r in core) == 161, 'CORE_ALL396_POSITIVE_COUNT_CHANGED')
    require(sum(r['native_fdr_pass_0_05'] == 'True' for r in core if r['original_interpretation_status'] == 'PRIMARY') == 153, 'CORE_PRIMARY_COUNT_CHANGED')
    require(sum(r['native_fdr_pass_0_05'] == 'True' for r in extension) == 603, 'EXTENSION_ALL1200_POSITIVE_COUNT_CHANGED')
    classes = rows['independent_validation_admissibility217_v4_1.tsv']
    require(len({r['pair_id'] for r in classes}) == 217, 'LOCKED217_UNIQUE_KEYS_REQUIRED')
    require(dict(collections.Counter(r['current_class'] for r in classes)) == {
        'QUALIFIED_EXTERNAL_OUTCOME_SIDE_VALIDATION': 23,
        'DIRECTIONALLY_CONCORDANT_BELOW_FROZEN_THRESHOLD': 18,
        'QC_INELIGIBLE_NOT_A_NULL_RG_TEST': 17,
        'NO_ELIGIBLE_EXTERNAL_SOURCE_NOT_A_NULL_RG_TEST': 159,
    }, 'ALL217_CLASSES_CHANGED')
    require(all(r['fully_independent_two_trait_replication'] == 'False' and r['shared_estimator_covariance_calibrated'] == 'False' and r['unestimated_or_QC_failed_is_biological_null'] == 'False' for r in classes), 'VALIDATION_SCOPE_QUALIFICATIONS_CHANGED')
    estimated = {r['pair_id'] for r in classes if r['current_native_estimate_present'] == 'True'}
    require(len(estimated) == 41 and estimated == {r['original_pair_id'] for r in rows['validation_native_full_precision_rg_v4.tsv']}, 'EXACT41_NATIVE_CLASSIFICATION_JOIN_REQUIRED')
    prior = rows['qualified_validation_prior_art23_v4_1.tsv']
    require({r['pair_id'] for r in prior} == {r['pair_id'] for r in classes if r['current_class'] == 'QUALIFIED_EXTERNAL_OUTCOME_SIDE_VALIDATION'}, 'EXACT23_PRIOR_ART_JOIN_REQUIRED')
    require(sorted(collections.Counter(r['prior_art_class'] for r in prior).values()) == [2, 6, 15] and all(r['new_primary_claim_admitted'] == 'False' for r in prior), 'PRIOR_ART_QUALIFICATIONS_CHANGED')
    for name, digest in observed.items():
        require(hashlib.sha256((PACKAGE / 'tables' / name).read_bytes()).hexdigest() == digest, 'TABLE_CHANGED_DURING_CHECK')
    print(json.dumps(dict(status='PORTABLE_AGGREGATE_TABLE_CHECK_PASS', table_sha256=observed,
        core_rg=396, extension_rg=1200, estimated_validation_rg=41, validation_candidates=217,
        sensitivity_estimates=62, fully_independent_two_trait_replications=0,
        new_primary_claims_admitted=0, SSD_or_GWAS_body_access=False,
        estimator_reexecution_or_native_vector_adjudication=False)))


if __name__ == '__main__':
    main()

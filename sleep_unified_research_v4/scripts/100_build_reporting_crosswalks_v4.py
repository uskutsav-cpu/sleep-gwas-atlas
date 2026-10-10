#!/usr/bin/env python3
"""Join sealed classification, native precision, cohort and prior-art evidence."""
import collections
import csv
import hashlib
import json
from pathlib import Path

P = Path(__file__).resolve().parents[1]
V1 = P.parent / 'sleep_unified_research_v1'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def rows(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def indexed(records, key):
    result = {r[key]: r for r in records}
    if len(result) != len(records):
        raise RuntimeError('DUPLICATE_REPORTING_JOIN_KEY: ' + key)
    return result

def write_table(path, records):
    with Path(path).open('x', newline='') as f:
        writer = csv.DictWriter(f, list(records[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)

def main():
    paths = dict(
        classes=P/'reviews/independent_whole_validation_classification217_v4.tsv',
        native=P/'tables/validation_native_full_precision_rg_v4.tsv',
        sources=V1/'REPLICATION_SOURCE_LEDGER.tsv',
        cohorts=V1/'reviews/cohort_independence_v1.tsv',
        prior=V1/'tables/historical_novelty_master_1200.tsv',
        prior_numeric=V1/'tables/historical_novelty_crosswalk.tsv',
        current_literature=P/'source_provenance/sleep_geneticist_literature_crosswalk_v4.tsv',
        native_completion=P/'logs/native190_root_independent_completion_addendum_v4.json')
    before = {str(path): sha(path) for path in [Path(__file__), *paths.values()]}
    classes = rows(paths['classes'])
    native = indexed(rows(paths['native']), 'original_pair_id')
    sources = indexed(rows(paths['sources']), 'pair_id')
    cohorts = indexed(rows(paths['cohorts']), 'pair_id')
    prior = indexed(rows(paths['prior']), 'pair_id')
    if len(classes) != 217 or len(native) != 41 or len(prior) != 1200:
        raise RuntimeError('LOCKED_FAMILY_COUNTS_REQUIRED')
    keys = {r['pair_id'] for r in classes}
    if len(keys) != 217 or keys != set(sources) or keys != set(cohorts):
        raise RuntimeError('EXACT_LOCKED217_JOIN_REQUIRED')
    expected = dict(QUALIFIED_EXTERNAL_OUTCOME_SIDE_VALIDATION=23,
        DIRECTIONALLY_CONCORDANT_BELOW_FROZEN_THRESHOLD=18,
        QC_INELIGIBLE_NOT_A_NULL_RG_TEST=17,
        NO_ELIGIBLE_EXTERNAL_SOURCE_NOT_A_NULL_RG_TEST=159)
    if dict(collections.Counter(r['independent_current_class'] for r in classes)) != expected:
        raise RuntimeError('UNCHANGED_COMPLETE217_CLASSIFICATION_REQUIRED')
    validation = []
    positives = []
    for r in classes:
        key = r['pair_id']; s = sources[key]; c = cohorts[key]; n = native.get(key)
        if r['independent_two_trait_replication'] != 'False' or r['all_classification_checks_pass'] != 'True' or (n is not None) != (r['native_estimate_present'] == 'True'):
            raise RuntimeError('CURRENT_INDEPENDENT_REVIEW_CLASSIFICATION_REQUIRED')
        if n and (n['rg'] != r['native_rg'] or n['p'] != r['native_p'] or n['independent_two_trait_replication'] != 'False'):
            raise RuntimeError('FULL_PRECISION_NATIVE_JOIN_CHANGED')
        record = dict(pair_id=key, sleep_trait=r['sleep_trait'],
            discovery_sleep_source_id=s['discovery_sleep_source_id'],
            discovery_sleep_cohort=c['discovery_sleep_cohort'],
            discovery_outcome_source=s['discovery_external_source'],
            external_validation_source_id=r['replication_source_id'],
            external_accession=s['replication_study_accession'],
            external_source_url=s['replication_source_url'],
            external_definition=s['replication_phenotype_definition'],
            phenotype_match=s['phenotype_match_status'],
            historical_class=r['historical_class'], current_class=r['independent_current_class'],
            current_native_estimate_present=r['native_estimate_present'],
            native_rg=n['rg'] if n else '', native_SE=n['se'] if n else '',
            native_P=n['p'] if n else '', CI_lower_95=n['ci_lower_95'] if n else '',
            CI_upper_95=n['ci_upper_95'] if n else '', exact_alpha=r['exact_alpha_0_05_over_217'],
            Dsleep_Vsleep_intersection=c['Dsleep_Vsleep_intersection'],
            Dsleep_Voutcome_intersection=c['Dsleep_Voutcome_intersection'],
            Doutcome_Vsleep_intersection=c['Doutcome_Vsleep_intersection'],
            Doutcome_Voutcome_intersection=c['Doutcome_Voutcome_intersection'],
            fully_independent_two_trait_replication='False',
            shared_estimator_covariance_calibrated='False',
            unestimated_or_QC_failed_is_biological_null='False',
            source_raw_preprocessing_replay='SEPARATE_PENDING_ORIGINAL13_REPLAY',
            classification_evidence=str(paths['classes']),
            native_estimate_evidence=str(paths['native']) if n else '',
            source_evidence=str(paths['sources']), cohort_evidence=str(paths['cohorts']))
        validation.append(record)
        if r['independent_current_class'] == 'QUALIFIED_EXTERNAL_OUTCOME_SIDE_VALIDATION':
            p = prior[key]
            positives.append(dict(pair_id=key, sleep_trait=r['sleep_trait'],
                extension_trait_id=p['extension_trait_id'], phenotype=p['phenotype_name'],
                native_external_validation_rg=n['rg'], native_external_validation_SE=n['se'],
                native_external_validation_P=n['p'], current_validation_class=r['independent_current_class'],
                prior_art_class=p['current_novelty_class'], comparability=p['current_pair_evidence_status'],
                prior_DOIs=p['prior_dois'], direct_comparator_sources=p['direct_comparator_sources'],
                exact_pair_search_scope=p['current_pair_specific_database_search'],
                literature_coverage=p['literature_coverage'],
                first_ever_relationship_established='False', new_primary_claim_admitted='False',
                unresolved_priority_is_novelty='False', remaining_review=p['remaining_review'],
                current_primary_study_crosswalk=str(paths['current_literature']),
                prior_pair_evidence=str(paths['prior']), native_evidence=str(paths['native'])))
    novelty_counts = dict(collections.Counter(r['prior_art_class'] for r in positives))
    if len(positives) != 23 or sorted(novelty_counts.values()) != [2,6,15]:
        raise RuntimeError('UNCHANGED_23_PAIR_NOVELTY_CLASSES_REQUIRED')
    if {path: sha(path) for path in before} != before:
        raise RuntimeError('REPORTING_JOIN_INPUT_CHANGED')
    outputs = [P/'tables/independent_validation_admissibility217_v4.tsv',
        P/'tables/qualified_validation_prior_art23_v4.tsv']
    write_table(outputs[0], validation); write_table(outputs[1], positives)
    if {path: sha(path) for path in before} != before:
        raise RuntimeError('REPORTING_INPUT_CHANGED_AFTER_PERSISTENCE')
    receipt = dict(status='EXACT_FAMILY_REPORTING_CROSSWALKS_COMPLETE',
        source_sha256=before, output_sha256={str(path):sha(path) for path in outputs},
        locked_candidates=217, native_estimates=41, classifications=expected,
        qualified_positive_prior_art_classes=novelty_counts,
        completed_fully_independent_two_trait_replications=0,
        first_ever_findings_established=0, new_primary_claims_admitted=0,
        source_replay_completion_claimed=False, scientific_assumptions_cleared=False)
    with (P/'logs/reporting_crosswalks_receipt_v4.json').open('x') as f:
        f.write(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(outputs=receipt['output_sha256'],novelty_counts=novelty_counts)))

if __name__ == '__main__':
    main()

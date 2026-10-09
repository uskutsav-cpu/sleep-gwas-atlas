"""Portable evidence invariants; not scientific estimator validation."""
import csv
import hashlib
import json
from pathlib import Path
import unittest

PACKAGE=Path(__file__).resolve().parents[1]

def load(path):return json.loads((PACKAGE/path).read_text())

def rows(path):
    with (PACKAGE/path).open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

class EvidenceContractsV3(unittest.TestCase):
    def test_source_admission_rejects_current_prostate(self):
        records=rows('tables/CURRENT_CORE_SOURCE_STATUS_V3.tsv')
        raw=[r for r in records if r['kind']=='core_raw'];archives=[r for r in records if r['kind']=='core_source_archive']
        self.assertEqual((len(raw),len(archives)),(45,43))
        self.assertEqual(sum(r['current_status']=='EXACT_HISTORICAL_HASH_AVAILABLE' for r in raw),45)
        self.assertEqual(sum(r['current_status']=='EXACT_HISTORICAL_HASH_AVAILABLE' for r in archives),42)
        r=next(r for r in archives if r['trait_id']=='prostate_cancer')
        self.assertNotEqual(r['expected_sha256'],r['public_candidate_sha256'])
        self.assertEqual(r['public_candidate_status'],'SHA_MISMATCH_PRESERVED')
        self.assertIn('UNRECOVERED',r['current_status'])

    def test_cumulative_missing_ledger(self):
        r=rows('tables/CUMULATIVE_MISSING_SOURCE_RECOVERY_V3.tsv')
        self.assertEqual(len(r),19)
        available=[x for x in r if x['current_status']!='UNRECOVERED']
        self.assertEqual(len(available),18)
        self.assertEqual(len({x['actual_sha256'] for x in available}),15)
        self.assertTrue(all(x['actual_sha256']==x['expected_sha256'] for x in available))
        self.assertEqual([x['trait_id'] for x in r if x['current_status']=='UNRECOVERED'],['prostate_cancer'])

    def test_MVP_complete_row_reconciliation(self):
        d=load('logs/mvp_preprocessing_worker_receipt_v3.json');e=load('logs/mvp_preprocessing_execution_receipt_v3.json')
        self.assertEqual(d['counts']['source_rows'],19703815);self.assertEqual(d['counts']['output_rows'],826027)
        self.assertEqual(sum(d['exclusive_exclusions'].values())+d['counts']['output_rows'],d['counts']['source_rows'])
        self.assertEqual(sum(d['retained_P_CI_relative_difference_buckets'].values()),d['counts']['output_rows'])
        self.assertEqual(sum(d['whole_file_log_CI_centering_SE_unit_buckets'].values()),d['counts']['source_rows'])
        self.assertTrue(e['input_hashes_unchanged']);self.assertTrue(e['execution_plan_hash_unchanged'])
        self.assertEqual(e['execution_plan_sha256'],sha(PACKAGE/'manifests/mvp_preprocessing_execution_plan_v3.json'))
        self.assertEqual(e['worker_receipt_sha256'],sha(PACKAGE/'logs/mvp_preprocessing_worker_receipt_v3.json'))
        self.assertFalse(d['h2_or_correlation_estimated']);self.assertEqual(d['pair_tests_admitted'],0)

    def test_independent_MVP_checks_do_not_admit_replication(self):
        d=load('reviews/mvp_preprocessing_stream_review_v3.json')
        self.assertEqual(d['status'],'INDEPENDENT_STREAMING_PREPROCESSING_PASS');self.assertEqual(d['check_count'],35)
        self.assertTrue(d['all_output_rows_streamed_and_CRC_verified']);self.assertTrue(d['raw_source_full_stream_CRC_verified'])
        self.assertEqual(d['pair_tests_admitted'],0);self.assertFalse(d['h2_or_rg_estimated'])

    def test_native_raw_content_equality_has_distinct_serialization_flag(self):
        r=rows('tables/NATIVE_CANCER_RAW_CONTENT_COMPARISON_V3.tsv')
        self.assertEqual({x['trait_id'] for x in r},{'breast_cancer','ovarian_cancer'})
        for x in r:
            self.assertEqual(x['old_decompressed_sha256'],x['new_decompressed_sha256'])
            self.assertEqual(x['gzip_bytes_identical'],'False')
            self.assertEqual(x['numeric_or_content_difference'],'NONE;EXACT_DECOMPRESSED_BYTES')
            self.assertEqual(x['full_GWAS_QC_or_estimator_reproduced'],'False')

    def test_HEAD_availability_is_not_body_reproduction(self):
        r=rows('tables/EXTENSION_LOCKED_VERSION_AVAILABILITY_V3.tsv')
        self.assertEqual(len(r),100);self.assertEqual(sum(int(x['expected_bytes']) for x in r),227610388647)
        self.assertTrue(all(x['status']=='EXACT_VERSION_HEAD_MATCH_BODY_NOT_REVERIFIED' and x['current_body_checksum_verified']=='False' for x in r))
        for x in r:self.assertEqual(sha(PACKAGE/x['header_proof']),x['header_sha256'])

    def test_figure_claims_bound_to_numerical_tables(self):
        d=load('manifests/mvp_QC_figure_receipt_v3.json')
        for p,s in {**d['tables'],**d['outputs']}.items():self.assertEqual(sha(PACKAGE/p),s)
        self.assertFalse(d['genetic_correlation_or_mechanism_claims'])

    def test_scientific_incompleteness_and_unmodified_guard(self):
        d=load('tables/CURRENT_EVIDENCE_STATUS_V3.json')
        for k in ['complete_core_396_native_reproduction','complete_extension_1200_native_reproduction','complete_validation_41_native_reproduction','HUMAN_MANUSCRIPT_AUTHORING_READY','scientific_goal_complete','shared_covariance_calibrated','new_primary_question_admitted']:
            self.assertFalse(d[k],k)
        self.assertEqual(d['both_trait_independent_successes'],0);self.assertEqual(d['new_pair_outcomes_estimated'],0)
        r=load('logs/native_resource_preflight_v3.json')
        self.assertEqual(r['internal_minimum_bytes'],3221225472);self.assertFalse(r['resource_gate_pass'])
        a=load('manifests/full_native_blocking_audit_v3.json')
        self.assertTrue(a['same_full_native_resource_condition_in_three_consecutive_goal_turns'])
        self.assertTrue(a['goal_not_complete'])

if __name__=='__main__':unittest.main()

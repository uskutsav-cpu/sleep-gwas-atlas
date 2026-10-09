import csv
import hashlib
import json
from pathlib import Path
import unittest

PACKAGE=Path(__file__).resolve().parents[1]
V1=PACKAGE.parent/'sleep_unified_research_v1'

def read(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))

def load(relative):
    return json.loads((PACKAGE/relative).read_text())

class RecoveryEvidence(unittest.TestCase):
    def test_missing_ledger_resolution_preserves_expected_sources(self):
        old=[r for r in read(V1/'tables/native_input_hash_checks.tsv') if r['status']=='MISSING']
        new=read(PACKAGE/'tables/SOURCE_RECOVERY_UPDATES.tsv')
        self.assertEqual(len(old),19)
        self.assertEqual({(r['kind'],r['trait_id'],r['path'],r['expected_sha256']) for r in old},
                         {(r['kind'],r['trait_id'],r['original_missing_path'],r['expected_sha256']) for r in new})
        resolved=[r for r in new if r['current_status']!='UNRECOVERED']
        self.assertEqual(len(resolved),16)
        for r in resolved:
            self.assertEqual(r['actual_sha256'],r['expected_sha256'])
            self.assertNotEqual(r['resolved_path'],r['original_missing_path'])
            self.assertEqual(r['original_missing_path_replaced'],'False')
        self.assertEqual({r['trait_id'] for r in new if r['current_status']=='UNRECOVERED'},
                         {'breast_cancer','prostate_cancer','ovarian_cancer'})

    def test_current_source_cardinality_and_proof(self):
        rows=read(PACKAGE/'tables/CURRENT_CORE_SOURCE_STATUS.tsv')
        raw=[r for r in rows if r['kind']=='core_raw']
        archive=[r for r in rows if r['kind']=='core_source_archive']
        self.assertEqual(len(raw),45)
        self.assertEqual(len({r['trait_id'] for r in raw}),45)
        self.assertTrue(all(r['current_sha256']==r['expected_sha256'] for r in raw))
        self.assertEqual(len(archive),43)
        self.assertEqual(sum(r['current_status']=='MISSING' for r in archive),3)
        for r in rows:
            self.assertTrue(r['source_release'] and r['phenotype_identity'] and r['effect_coding'])

    def test_ranged_source_admission_and_independent_interval_coverage(self):
        aggregate=load('logs/exact_lipid_ranged_acquisition_receipt_v2.json')
        self.assertTrue(aggregate['all_exact_expected_hashes'])
        self.assertEqual({r['trait_id'] for r in aggregate['sources']},{'hdl','ldl','triglycerides'})
        self.assertEqual(aggregate['total_source_bytes'],6844892917)
        for source in aggregate['sources']:
            self.assertEqual(source['actual_bytes'],source['expected_bytes'])
            self.assertEqual(source['actual_sha256'],source['expected_sha256'])
            detail=load('logs/'+source['trait_id']+'_ranged_acquisition_receipt_v2.json')
            cursor=0
            for chunk in sorted(detail['chunks'],key=lambda r:r['start']):
                self.assertEqual(chunk['start'],cursor)
                self.assertEqual(chunk['actual_bytes'],chunk['end']-chunk['start']+1)
                self.assertEqual(chunk['status'],'PASS_EXACT_RANGE')
                self.assertEqual(chunk['headers']['content-range'],f"bytes {chunk['start']}-{chunk['end']}/{source['expected_bytes']}")
                cursor=chunk['end']+1
            self.assertEqual(cursor,source['expected_bytes'])
        failed=load('logs/hdl_source_acquisition_receipt_v2.json')
        self.assertEqual(failed['status'],'FAILED_PARTIAL_OR_MISMATCH_PRESERVED')
        self.assertIsNone(failed['actual_sha256'])

    def test_original_native_prefilter_stage_exactness(self):
        historical={r['trait_id']:r for r in read(V1/'tables/native_input_hash_checks.tsv') if r['kind']=='core_prefiltered'}
        counts={'hdl':(46150908,1223518),'ldl':(47006483,1223471),'triglycerides':(47196261,1223515)}
        for trait,(total,retained) in counts.items():
            r=load('logs/'+trait+'_native_prefilter_receipt_v2.json')
            self.assertEqual(r['status'],'NATIVE_PREFILTER_EXACT_HISTORICAL_BYTES_AND_COUNTS')
            self.assertEqual(r['comparison']['source_rows'],[total,total])
            self.assertEqual(r['comparison']['retained_rows'],[retained,retained])
            self.assertEqual(r['output_sha256'],historical[trait]['expected_sha256'])
            self.assertEqual(r['inputs_sha256_before'],r['inputs_sha256_after'])
            self.assertFalse(r['complete_GWAS_QC_or_estimator_reproduction'])
            plan_path=PACKAGE/'manifests'/(trait+'_prefilter_execution_plan_v2.json')
            self.assertEqual(hashlib.sha256(plan_path.read_bytes()).hexdigest(),r['resource_plan_sha256'])

    def test_prospective_definition_without_new_outcomes(self):
        r=load('manifests/protocol_amendment_v3_receipt.json')
        self.assertTrue(r['public_clinical_definition_gate_resolved'])
        self.assertFalse(r['source_h2_or_pair_outcomes_examined'])
        self.assertEqual(r['new_pair_tests_admitted'],0)
        for path,expected in r['bound_files'].items():
            self.assertEqual(hashlib.sha256((PACKAGE.parent/path).read_bytes()).hexdigest(),expected)
        definition=load('source_definition_review/provenance_v2_definition_contract.json')
        self.assertEqual(definition['n_cases']+definition['n_controls'],408138)
        self.assertEqual(sum(len(c) for c in definition['code_lists'].values()),18)
        self.assertFalse(definition['control_rule']['related_phecode_exclusions_applied'])

    def test_scientific_completion_is_distinct_from_integrity(self):
        r=load('tables/CURRENT_EVIDENCE_STATUS_V2.json')
        self.assertEqual(r['raw_sources_exact_available'],45)
        self.assertEqual(r['source_containers_exact_available'],40)
        self.assertEqual(r['native_prefilters_exact'],3)
        self.assertEqual(r['both_trait_independent_successes'],0)
        self.assertEqual(r['new_clinical_pair_outcomes_examined'],0)
        self.assertFalse(r['HUMAN_MANUSCRIPT_AUTHORING_READY'])
        self.assertFalse(r['complete_extension_1200_native_reproduction'])
        self.assertEqual(r['full_native_resource_preflight']['internal_minimum_bytes'],3*1024**3)
        self.assertFalse(r['full_native_resource_preflight']['full_native_job_launched'])

if __name__=='__main__':
    unittest.main()

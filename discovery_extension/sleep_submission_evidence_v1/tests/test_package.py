import importlib.util
import unittest
from pathlib import Path
P=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('validate_package',P/'scripts/validate_package.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)

class EvidenceContractTests(unittest.TestCase):
    def test_duplicate_identity_cannot_shrink_universe(self):
        with self.assertRaises(ValueError):v.indexed([{'id':'a'},{'id':'a'}],'id')
    def test_discovery_and_replication_evidence_remain_consistent(self):
        self.assertEqual(v.check()['status'],'PASS_PACKAGE_CONSISTENCY')
    def test_missing_replication_effects_are_not_zero(self):
        rows=v.read(P/'tables/replication_family_217.tsv')
        unavailable=[r for r in rows if r['replication_class']=='NO_INDEPENDENT_DATASET']
        self.assertEqual(len(unavailable),159)
        self.assertTrue(all(r['replication_rg']=='NA' and r['replication_se']=='NA' for r in unavailable))
    def test_no_two_trait_independence_overclaim(self):
        rows=v.read(P/'tables/replicated_23.tsv')
        self.assertTrue(all(r['audit_replication_class']=='EXTERNAL_OUTCOME_SIDE_REPLICATION' for r in rows))
    def test_novelty_priority_is_not_inferred_from_no_match(self):
        rows=v.read(P/'tables/NOVELTY_MASTER_1200.tsv')
        self.assertTrue(all(r['adequately_supported_novel_result'] not in ('True','TRUE','YES') for r in rows))

if __name__=='__main__':unittest.main()

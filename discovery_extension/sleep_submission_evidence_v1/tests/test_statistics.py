"""Numerical contracts and failure-mode checks; source-free, no scientific replay."""
import importlib.util
import math
import unittest
from pathlib import Path

BASE=Path(__file__).resolve().parents[1]

def load(name):
    spec=importlib.util.spec_from_file_location(name,BASE/'scripts'/f'{name}.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

v=load('verify_statistics');ind=load('independent_statistics')

class StatisticsTests(unittest.TestCase):
    def test_known_normal_probabilities(self):
        for z,p in [(0,1),(1,0.31731050786291415),(1.959963984540054,.05),(3,0.0026997960632601913)]:
            self.assertAlmostEqual(v.normal_p(z),p,places=14)
            self.assertAlmostEqual(ind.tail(z),p,places=14)
    def test_independent_extreme_tail_agreement(self):
        for z in (.01,.1,1,2,5,8,12,15,20,30):
            self.assertTrue(math.isclose(ind.tail(z),v.normal_p(z),rel_tol=2e-12,abs_tol=1e-300))
    def test_signed_tails_symmetric(self):
        self.assertEqual(v.normal_p(-5),v.normal_p(5))
        self.assertEqual(ind.tail(-5),ind.tail(5))
    def test_bh_unsorted_ties_zero_one(self):
        p=[.2,.01,.01,0,1]
        expected=[.25,1/60,1/60,0,1]
        for method in (v.bh,ind.adjustment):
            for x,y in zip(method(p),expected):self.assertAlmostEqual(x,y)
    def test_bh_rejects_incomplete_or_invalid_family(self):
        for p in ([],[float('nan')],[float('inf')],[-.1],[1.1]):
            for method in (v.bh,ind.adjustment):
                with self.assertRaises(ValueError):method(p)
    def test_bh_not_substituted_by_selected_subset(self):
        self.assertNotEqual(v.bh([.03,.8,.9])[0],v.bh([.03])[0])
    def test_nonfinite_tail_rejected(self):
        for z in (float('nan'),float('inf')):
            for method in (v.normal_p,ind.tail):
                with self.assertRaises(ValueError):method(z)
    def test_missing_preserved_as_missing(self):
        for x in ('NA','',None,'NOT_AVAILABLE'):self.assertIsNone(v.numeric(x))
        self.assertEqual(v.numeric('0'),0)
    def test_duplicate_identity_rejected(self):
        with self.assertRaises(ValueError):v.unique([{'id':'x'},{'id':'x'}],'id')
    def test_invalid_standard_error_rejected(self):
        for se in (0,-1,float('nan')):
            with self.assertRaises(ValueError):v.estimates(.1,se)
    def test_unconstrained_ci_not_truncated(self):
        _,_,lo,hi=v.estimates(.95,.1)
        self.assertGreater(hi,1);self.assertLess(lo,.95)
    def test_covariance_changes_difference_uncertainty(self):
        self.assertGreater(v.difference(.5,.1,.2,.1,-.5)[1],v.difference(.5,.1,.2,.1,.5)[1])
        self.assertLess(v.difference(.5,.1,.2,.1,.5)[4],v.difference(.5,.1,.2,.1,0)[4])
    def test_difference_direction_and_null(self):
        self.assertEqual(v.difference(.5,.1,.5,.2)[4],1)
        self.assertLess(v.difference(.5,.1,.2,.1)[0],0)
    def test_impossible_covariance_rejected(self):
        for rho in (-1.01,1.01):
            with self.assertRaises(ValueError):v.difference(.5,.1,.2,.1,rho)
        with self.assertRaises(ValueError):v.difference(.5,.1,.2,.1,1)
    def test_generated_universes_unique_complete(self):
        root=BASE.parents[1]
        g=v.read(BASE/'tables/global_1200.tsv');r=v.read(BASE/'tables/replication_family_217.tsv')
        p=v.read(root/'discovery_extension/config/candidate_traits.tsv')
        s=[x for x in v.read(root/'config/analysis_panel.tsv') if x['domain']=='sleep']
        self.assertEqual(v.unique(g,'pair_id'),{x['trait_id']+'__'+y['extension_trait_id'] for x in s for y in p})
        self.assertEqual(v.unique(r,'pair_id'),v.unique(v.read(root/'discovery_extension/config/replication_manifest.tsv'),'pair_id'))
    def test_no_missing_replicate_converted_to_zero(self):
        for r in v.read(BASE/'tables/replication_family_217.tsv'):
            if r['replication_class'] in ('UNDERPOWERED','NO_INDEPENDENT_DATASET'):
                self.assertEqual(r['replication_rg'],'NA');self.assertEqual(r['replication_p'],'NA')
    def test_independence_qualification_not_full_replication(self):
        for r in v.read(BASE/'tables/replicated_23.tsv'):
            self.assertEqual(r['sleep_GWAS_reused'],'True')
            self.assertEqual(r['audit_replication_class'],'EXTERNAL_OUTCOME_SIDE_REPLICATION')
    def test_actigraphy_measurement_classification_matches_definition(self):
        root=BASE.parents[1]
        core={x['trait_id']:x for x in v.read(root/'config/analysis_panel.tsv') if x['domain']=='sleep'}
        for r in v.read(BASE/'tables/global_1200.tsv'):
            if 'actigraphy' in core[r['sleep_trait']]['phenotype_definition'].lower():
                self.assertEqual(r['sleep_measurement'],'ACCELEROMETER')
    def test_negative_secondary_queue_is_exact_discovery_complement(self):
        root=BASE.parents[1]
        nulls={r['pair_id'] for r in v.read(BASE/'tables/global_1200.tsv') if r['historical_significant']=='False'}
        queue=v.read(root/'discovery_extension/results/local/local_analysis_queue.tsv')
        self.assertEqual(nulls,{r['pair_id'] for r in queue if r['selection_status']=='PENDING_GLOBAL_NULL_CURATION'})
    def test_correction_does_not_shrink_to_available_tests(self):
        rows=v.read(BASE/'tables/replication_family_217.tsv')
        for r in rows:self.assertEqual(float(r['replication_alpha_exact']),.05/len(rows))

if __name__=='__main__':unittest.main()

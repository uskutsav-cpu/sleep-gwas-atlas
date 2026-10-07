"""Synthetic numerical tests, deliberately separate from empirical results."""
import importlib.util
from pathlib import Path
import unittest

p=Path(__file__).resolve().parents[1]/'scripts/audit_evidence.py'
s=importlib.util.spec_from_file_location('continuation_audit',p)
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

class NumericalChecks(unittest.TestCase):
    def test_bh_unsorted_ties_and_full_denominator(self):
        self.assertEqual(m.bh([.9,.01,.02,.02]),[.9,.02666666666666667,.02666666666666667,.02666666666666667])
    def test_bh_against_independent_scipy(self):
        from scipy.stats import false_discovery_control
        ps=[.004,.97,.03,.003,.03,.11,.2]
        for a,b in zip(m.bh(ps),false_discovery_control(ps)):self.assertAlmostEqual(a,b,places=14)
    def test_prior_reweight_matches_odds_scaling(self):
        v=[.01,.02,.03,.04,.9];r=m.reweight(v,.1)
        self.assertAlmostEqual(r[4],.09/.19)
        self.assertAlmostEqual(sum(r),1)
        with self.assertRaises(ValueError):m.reweight([0,0,0,0,-1],.1)
    def test_overlap_is_pair_specific_inclusive(self):
        c={'candidate_locus_id':'placo_insomnia__adhd_chr5_10_20','pair_id':'insomnia__adhd'}
        b={'pair_id':'insomnia__adhd','chr':'5','start':'20','end':'30'}
        self.assertTrue(m.overlap(c,b));b['pair_id']='insomnia__mdd';self.assertFalse(m.overlap(c,b))

if __name__=='__main__':unittest.main()

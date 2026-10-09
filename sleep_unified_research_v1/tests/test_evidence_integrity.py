"""Evidence invariants. Passing these tests does not complete scientific validation."""
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import unittest

PACKAGE=Path(__file__).resolve().parents[1]
def rows(relative):
    with (PACKAGE/relative).open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def load(relative):return json.loads((PACKAGE/relative).read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bh(values):
    order=sorted(range(len(values)),key=values.__getitem__);out=[None]*len(values);minimum=1.
    for rank in range(len(order),0,-1):
        index=order[rank-1];minimum=min(minimum,values[index]*len(order)/rank);out[index]=minimum
    return out

class EvidenceIntegrity(unittest.TestCase):
    def test_all_18_checkpoint_artifacts_are_exact(self):
        recovered=rows('ORIGINAL_CORE_RECOVERY.tsv');self.assertEqual(len(recovered),18)
        for r in recovered:
            self.assertEqual(r['status'],'EXACT_SHA256_RECOVERED')
            self.assertEqual(sha(PACKAGE/r['recovered_path']),r['expected_sha256'])
    def test_separate_global_families_and_frozen_BH(self):
        for name,size,key,q,total,primary in [('tables/original_core_396.tsv',396,'disease_trait','fdr',161,153),('tables/original_extension_1200.tsv',1200,'extension_trait_id','extension_fdr',603,603)]:
            data=rows(name);self.assertEqual(len(data),size)
            self.assertEqual(len({(r['sleep_trait'],r[key]) for r in data}),size)
            calculated=bh([float(r['p']) for r in data])
            self.assertEqual(sum(v<.05 for v in calculated),total)
            for r,v in zip(data,calculated):self.assertEqual(float(r[q]),v)
            self.assertEqual(sum(v<.05 and r.get('analysis_tier','PRIMARY_PHASE1')=='PRIMARY_PHASE1' for r,v in zip(data,calculated)),primary)
    def test_replication_family_independence_and_QC(self):
        data=rows('REPLICATION_RESULTS.tsv');self.assertEqual(len(data),217)
        self.assertEqual(len({r['pair_id'] for r in data}),217)
        self.assertEqual(sum(r['historical_replication_class']=='REPLICATED' for r in data),23)
        self.assertEqual(sum(r['current_QC_reason']=='H2_Z_LT4' for r in data),9)
        self.assertEqual(sum(r['current_QC_reason']=='H2_INTERCEPT_GT1P2' for r in data),8)
        self.assertEqual(sum(r['current_participant_overlap_evidence'].startswith('COHORT_SOURCE_DISTINCT') for r in data),58)
        self.assertTrue(all(r['both_trait_independent_replication']=='false' for r in data))
        estimated=[r for r in data if r['replication_rg'] not in ['','NA']]
        self.assertEqual(len(estimated),41)
        self.assertEqual(sum(float(r['replication_p'])<.05/217 and r['direction_concordant']=='True' for r in estimated),23)
        self.assertEqual(sum(float(r['heterogeneity_p'])<.05 for r in estimated),16)
        self.assertEqual(sum(float(r['heterogeneity_p'])<.05 and r['historical_replication_class']=='REPLICATED' for r in estimated),7)
    def test_native_pilot_ratio_jackknife(self):
        p=PACKAGE/'native/core_pilot_v1';e=load('native/core_pilot_v1/rg_insomnia__bmi.full_precision.json')['estimates'][0]
        delete={}
        for name in ['hsq1','hsq2','gencov']:
            files=list(p.glob('*.'+name+'.delete'));self.assertEqual(len(files),1)
            delete[name]=[float(x) for x in files[0].read_text().splitlines()];self.assertEqual(len(delete[name]),200)
        rg=e['gencov']['tot']/math.sqrt(e['hsq1']['tot']*e['hsq2']['tot'])
        dr=[c/math.sqrt(a*b) for a,b,c in zip(delete['hsq1'],delete['hsq2'],delete['gencov'])]
        mean=math.fsum(dr)/200
        se=math.sqrt(199/200*math.fsum((v-mean)**2 for v in dr))
        self.assertTrue(math.isclose(rg,e['rg_ratio'],rel_tol=1e-13))
        self.assertTrue(math.isclose(se,e['rg_se'],rel_tol=1e-13))
        self.assertTrue(math.isclose(math.erfc(abs(rg/se)/math.sqrt(2)),e['p'],rel_tol=1e-11))
    def test_independent_review_receipts(self):
        r=load('reviews/independent_numerical_v1.json');self.assertTrue(r['all_arithmetic_checks_pass'])
        self.assertEqual(r['independent_code_sha256'],sha(PACKAGE/'reviews/independent_numerical_v1.py'))
        r=load('reviews/extension_log_concordance_v1.json')
        self.assertEqual(r['status'],'PASS');self.assertFalse(r['mismatches']);self.assertEqual(r['original_rg_pairs_checked'],1200)
        self.assertEqual(r['original_h2_traits_checked'],100);self.assertEqual(r['frozen_1200_bh_significant_count'],603)
        r=load('reviews/cohort_independence_v1.json');self.assertEqual(r['historical_independent_both_trait_replications'],0)
    def test_archive_recovery_is_not_native_replay(self):
        r=load('logs/archived_extension_recovery_receipt_v1.json')
        self.assertEqual(r['recovered_receipt_pinned_artifacts'],227);self.assertFalse(r['missing_receipt_pinned_artifacts'])
        self.assertEqual(r['extension_munged_sources'],100);self.assertEqual(r['replication_munged_sources'],13)
        r=load('logs/progress_evidence_receipt_v1.json');self.assertFalse(r['full_native_scientific_reproduction_complete'])
        self.assertEqual((r['native_core_pair_count'],r['native_extension_pair_count']),(1,0))
    def test_native_runner_rejects_nonfinite_success(self):
        p=PACKAGE/'scripts/04_native_reproduction_runner.py';s=importlib.util.spec_from_file_location('runner',p)
        module=importlib.util.module_from_spec(s);s.loader.exec_module(module)
        good={'status':'NATIVE_ESTIMATE_RETURNED','rg_ratio':.2,'rg_se':.1,'z':2.,'p':.0455}
        for key in ['hsq1','hsq2','gencov']:good[key]={'tot':.1,'tot_se':.01,'intercept':1.,'intercept_se':.02}
        self.assertTrue(module.finite_estimates('rg',[good],1))
        self.assertFalse(module.finite_estimates('rg',[{**good,'rg_ratio':math.nan}],1))
        self.assertFalse(module.finite_estimates('rg',[{**good,'p':None}],1))
        self.assertFalse(module.finite_estimates('rg',[good],2))
    def test_protocol_identity_and_prospective_stop(self):
        v1=sha(PACKAGE/'FROZEN_NEW_ANALYSIS_PROTOCOL.md')
        self.assertEqual(load('manifests/native_reproduction_jobs_v1.json')['protocol_sha256'],v1)
        self.assertEqual(load('logs/mvp_insomnia_acquisition_receipt_v1.json')['protocol_sha256'],v1)
        a=load('manifests/protocol_amendment_v2_receipt.json')
        self.assertEqual(a['v1_sha256'],v1);self.assertEqual(a['v2_sha256'],sha(PACKAGE/'FROZEN_NEW_ANALYSIS_PROTOCOL_v2.md'))
        self.assertFalse(a['new_pair_tests_admitted']);self.assertFalse(a['new_pair_rg_outcomes_accessed'])
    def test_native_plan_retains_cardinality(self):
        plan=load('manifests/native_reproduction_jobs_v1.json');self.assertEqual(len(plan['jobs']),190)
        for stage,total in [('core',396),('extension',1200),('validation',41)]:
            self.assertEqual(sum(j['estimates'] for j in plan['jobs'] if j['stage']==stage and j['kind']=='rg'),total)
    def test_figures_bind_to_actual_tables(self):
        m=load('manifests/figure_manifest_v1.json');self.assertEqual(len(m['figures']),7)
        for p,h in {**m['input_sha256'],**m['output_sha256']}.items():self.assertEqual(sha(PACKAGE/p),h)
        self.assertFalse(m['new_biological_claims']);self.assertFalse(m['mechanism_figure_admitted'])

if __name__=='__main__':unittest.main()

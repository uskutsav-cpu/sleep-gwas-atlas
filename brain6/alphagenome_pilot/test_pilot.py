#!/usr/bin/env python3
"""Independent consistency checks against the frozen inputs and cached references."""
import collections,csv,gzip,hashlib,importlib.util,json,math,os,subprocess,sys,tempfile,types,unittest
from unittest import mock
from pathlib import Path
from pyliftover import LiftOver

D=Path(__file__).resolve().parent; ROOT=D.parents[1]; OUT=D/'results'; UP=ROOT/'brain6/results/brain6_exploratory_finemap_coloc_v2'; CHAIN=Path('/Volumes/Extreme SSD/brain6-work/independent-ld-reference-v1/hg19ToHg38.over.chain.gz')
def rows(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
class PilotTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.v=rows(OUT/'variant_evidence.tsv');cls.r=rows(OUT/'regions.tsv');cls.p=json.loads((OUT/'provenance.json').read_text())
    def test_binding_and_integrity(self):
        for name,expected in self.p['output_sha256'].items(): self.assertEqual(sha(OUT/name),expected)
        for name,expected in self.p['sequence_sha256'].items(): self.assertEqual(sha(OUT/'sequence_cache'/name),expected)
        self.assertEqual(sha(CHAIN),self.p['chain_sha256'])
        self.assertEqual(self.p['n_regions'],6);self.assertEqual(len(self.v),425)
    def test_region_and_membership_completeness(self):
        ids={r['candidate_locus_id'] for r in self.r};self.assertEqual(len(ids),6)
        upstream=rows(UP/'fine_mapping_credible_sets.tsv')
        expected={(r['candidate_locus_id'],r['snp'],r['trait_id']) for r in upstream if r['candidate_locus_id'] in ids}
        got={(r['candidate_locus_id'],r['snp'],trait) for r in self.v for trait in ['insomnia','adhd'] if r[trait+'_credible_set_member']=='YES'}
        self.assertEqual(expected,got)
        self.assertEqual(len({(r['candidate_locus_id'],r['snp']) for r in self.v}),len(self.v))
        counts=collections.Counter(r['candidate_locus_id'] for r in self.v)
        for r in self.r:self.assertEqual(counts[r['candidate_locus_id']],int(r['n_union_credible_variants']))
    def test_each_liftover_and_reference_base(self):
        lift=LiftOver(str(CHAIN)); comp=str.maketrans('ACGT','TGCA')
        for r in self.v:
            self.assertEqual(r['alphagenome_status'],'NOT_RUN_NO_API_KEY')
            if r['mapping_status']!='READY_UNIQUE_REFERENCE_VERIFIED':continue
            chrom=r['chromosome'];p=int(r['pos_grch37_1based'])-1
            hits=lift.convert_coordinate(chrom,p);self.assertEqual(len(hits),1)
            c,p38,strand,_=hits[0];self.assertEqual(c,chrom);self.assertEqual(p38+1,int(r['pos_grch38_1based']))
            def cached(build,pos):
                start=pos//100000*100000; obj=json.loads((OUT/'sequence_cache'/f'{build}_{chrom}_{start}_{start+100000}.json').read_text());return obj['dna'][pos-start].upper()
            ref37=cached('hg19',p);ref38=cached('hg38',p38)
            self.assertEqual(ref37,r['ref_grch37']);self.assertEqual(ref38,r['ref_grch38'])
            self.assertEqual(ref38,ref37.translate(comp) if strand=='-' else ref37)
            self.assertNotEqual(ref38,r['alt_grch38'])
            self.assertEqual(r['atlas_variant'],f"{chrom}:{p38+1}:{ref38}>{r['alt_grch38']}")
    def test_scoring_flatten_preserves_metadata_and_null(self):
        spec=importlib.util.spec_from_file_location('atlas_score',D/'atlas_score.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
        import numpy as np
        class Frame:
            def __init__(self,records):self.records=records;self.iloc=self
            def __getitem__(self,i):return self.records[i]
        class A:
            n_obs=1;n_vars=2;obs=Frame([{'gene_id':'ENSG1'}]);var=Frame([{'name':'brain1'},{'name':'brain2'}]);X=np.array([[.4,float('nan')]]);layers={'quantiles':np.array([[.7,float('nan')]])}
        out=m.flatten({'RNA_SEQ':A()});self.assertEqual(len(out),2);self.assertEqual(out[0]['observation']['gene_id'],'ENSG1');self.assertEqual(out[0]['track']['name'],'brain1');self.assertAlmostEqual(out[0]['quantile_score'],.7);self.assertIsNone(out[1]['raw_score'])
    def test_join_preserves_scored_and_unqueried_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp); run=base/'run'; (run/'scores').mkdir(parents=True)
            e=base/'evidence.tsv'; e.write_text('candidate_locus_id\tatlas_variant\nfirst\tchr5:1:A>C\nsecond\tchr5:2:G>T\n')
            (run/'run_config.json').write_text(json.dumps({'evidence_sha256':sha(e)}))
            score={'scorer':'RNA_SEQ','observation':{'gene_id':'ENSG1'},'track':{'name':'brain'},'raw_score':0.4,'quantile_score':None}
            with gzip.open(run/'scores'/'one.json.gz','wt') as f:json.dump({'variant':'chr5:1:A>C','scores':[score]},f)
            (run/'query_status.tsv').write_text('variant\tstatus\tn_score_rows\tfile\terror_type\nchr5:1:A>C\tSCORED\t1\tscores/one.json.gz\t\n')
            subprocess.run([sys.executable,str(D/'integrate.py'),'--evidence',str(e),'--run',str(run),'--out',str(base/'joined')],check=True,capture_output=True)
            got=rows(base/'joined'/'variant_molecular_evidence.tsv')
            self.assertEqual([x['atlas_query_status'] for x in got],['SCORED','NOT_QUERIED'])
            self.assertEqual(got[0]['raw_score'],'0.4');self.assertEqual(got[1]['raw_score'],'')
    def test_avi_join_and_independent_chr5_check(self):
        directory=OUT/'avi_mirror'
        if not directory.exists(): self.skipTest('AVI mirror has not been queried')
        receipt=json.loads((directory/'avi_provenance.json').read_text())
        result=rows(directory/'variant_evidence_avi.tsv')
        self.assertEqual(sha(directory/'variant_evidence_avi.tsv'),receipt['output_sha256'])
        self.assertEqual(sha(OUT/'variant_evidence.tsv'),receipt['evidence_sha256'])
        self.assertEqual(len(result),len(self.v));self.assertEqual(len(result),425)
        for original,joined in zip(self.v,result):
            self.assertTrue(all(joined[k]==v for k,v in original.items()))
            self.assertEqual(joined['avi_status'],'SCORED_MIRROR_EXPLORATORY')
            self.assertLessEqual(float(joined['avi_phred_lower']),float(joined['avi_phred_upper']))
        target=next(x for x in result if x['snp']=='rs2431108' and 'chr5_103447968' in x['candidate_locus_id'])
        self.assertEqual(target['atlas_variant'],'chr5:104612267:T>C')
        # Broad's independent lookup displays raw=0.559, PHRED=20.00 rounded.
        self.assertAlmostEqual(float(target['avi_raw_score_decimal']),0.559,delta=.0005)
        self.assertAlmostEqual(float(target['avi_phred_lower']),20.00,delta=.005)
    def test_descriptive_rank_comparison_and_gates(self):
        directory=OUT/'avi_mirror'; receipt=json.loads((directory/'ranking_provenance.json').read_text())
        for name,digest in receipt['outputs_sha256'].items():self.assertEqual(sha(directory/name),digest)
        comparison=rows(directory/'ranking_overlap.tsv');ranks=rows(directory/'variant_ranks.tsv')
        self.assertEqual(len(comparison),12);self.assertEqual(len(ranks),425)
        for r in comparison:
            if r['trait_gate']!='EXPLORATORY_ABF_CREDIBLE_SET':
                self.assertEqual(r['analysis_label'],'NOT_ESTIMATED_TRAIT_GATE')
                self.assertEqual(r['n_overlap'],'')
        target=next(x for x in ranks if x['snp']=='rs2431108' and 'chr5_103447968' in x['candidate_locus_id'])
        self.assertEqual((target['avi_rank_within_union'],target['insomnia_rank_within_union'],target['adhd_rank_within_union']),('1','1','3'))
    def test_broad_crosscheck_bound_to_each_region(self):
        directory=OUT/'avi_mirror'; path=directory/'broad_crosscheck_provenance.json'
        if not path.exists():self.skipTest('Cross-check is not complete')
        receipt=json.loads(path.read_text());checks=rows(directory/'broad_crosscheck.tsv')
        self.assertEqual(sha(directory/'variant_evidence_avi.tsv'),receipt['source_evidence_sha256'])
        self.assertEqual(sha(directory/'broad_crosscheck.tsv'),receipt['output_sha256'])
        self.assertEqual(len(checks),6);self.assertEqual(len({r['candidate_locus_id'] for r in checks}),6)
        evidence=rows(directory/'variant_evidence_avi.tsv')
        for check in checks:
            candidates=[x for x in evidence if x['candidate_locus_id']==check['candidate_locus_id']]
            expected=sorted(candidates,key=lambda x:(-int(x['avi_raw_e5']),x['atlas_variant']))[0]
            self.assertEqual(check['atlas_variant'],expected['atlas_variant'])
            if check['status']=='AGREES_TO_ROUNDING':
                self.assertLessEqual(abs(float(check['broad_raw'])-float(check['mirror_raw'])),.000501)
                self.assertLessEqual(float(check['mirror_phred_low'])-.00501,float(check['broad_phred']))
                self.assertGreaterEqual(float(check['mirror_phred_high'])+.00501,float(check['broad_phred']))
    def test_measured_eqtl_context_and_missing_states(self):
        directory=OUT/'avi_mirror'; meta=json.loads((directory/'molecular_context_provenance.json').read_text())
        for name,digest in meta['outputs_sha256'].items():self.assertEqual(sha(directory/name),digest)
        qtl=rows(directory/'lead_eqtl_context.tsv');status=rows(directory/'lead_eqtl_query_status.tsv')
        self.assertEqual(len(qtl),22);self.assertEqual(len(status),12)
        case=[r for r in qtl if r['lead_snp']=='rs2431108' and 'chr5_103447968' in r['locus_id']]
        self.assertEqual(len(case),2);self.assertEqual({r['nominal_pvalue'] for r in case},{'0.291197','0.0989162'})
        self.assertTrue(all(r['avi_match_status']=='EXACT_LOCUS_RSID_POSITION_UNORDERED_ALLELES' for r in case))
        self.assertTrue(any(r['status']=='NO_EXACT_LEAD_QTL_ROW' for r in status))
        self.assertTrue(any(r['avi_match_status']=='LEAD_NOT_IN_SCORED_CREDIBLE_SET_UNION' for r in status))
    def test_mocked_atlas_workflow_tracks_partial_results_and_secrets(self):
        import numpy as np
        import pandas as pd
        spec=importlib.util.spec_from_file_location('pilot_atlas',D/'atlas_score.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
        class Data:
            n_obs=1;n_vars=1
            obs=pd.DataFrame([{'gene_id':'ENSG1'}]);var=pd.DataFrame([{'name':'brain cortex'}])
            X=np.array([[0.42]]);layers={'quantiles':np.array([[0.7]])}
        class Client:
            def scorer_metadata(self):
                return {'RNA_SEQ':types.SimpleNamespace(is_signed=True,track_metadata=pd.DataFrame([{'ontology_curie':'UBERON:0001870'}]))}
            def query_variant(self,v,requested_scorers,ontology_terms):
                self_args=(requested_scorers,ontology_terms)
                if self_args!= (['RNA_SEQ'],['UBERON:0001870']):raise AssertionError('Bad query filters')
                if v=='chr5:2:G>T':raise RuntimeError('FAKE_SECRET must not be recorded')
                return {'RNA_SEQ':Data()}
        fake_atlas=types.ModuleType('alphagenome.atlas');fake_atlas.atlas=types.SimpleNamespace(create=lambda key,timeout:Client())
        fake_data=types.ModuleType('alphagenome.data');fake_data.genome=types.SimpleNamespace(Variant=types.SimpleNamespace(from_str=lambda s:s))
        modules={'alphagenome':types.ModuleType('alphagenome'),'alphagenome.atlas':fake_atlas,'alphagenome.data':fake_data}
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);e=p/'evidence.tsv';o=p/'run'
            e.write_text('candidate_locus_id\tmapping_status\tatlas_variant\nA\tREADY_UNIQUE_REFERENCE_VERIFIED\tchr5:1:A>C\nA\tREADY_UNIQUE_REFERENCE_VERIFIED\tchr5:2:G>T\nA\tHOLD_PALINDROMIC_ALLELES\t\n')
            argv=['atlas_score.py','score','--evidence',str(e),'--out',str(o),'--scorer','RNA_SEQ','--ontology','UBERON:0001870','--region','all']
            with mock.patch.dict(sys.modules,modules),mock.patch.dict(os.environ,{'ALPHAGENOME_API_KEY':'FAKE_SECRET'}),mock.patch.object(sys,'argv',argv),mock.patch.object(mod.importlib.metadata,'version',return_value='0.10.0'):
                mod.main()
            status=rows(o/'query_status.tsv');self.assertEqual([r['status'] for r in status],['SCORED','QUERY_FAILED'])
            with gzip.open(o/status[0]['file'],'rt') as f: result=json.load(f)
            self.assertEqual(result['scores'][0]['observation']['gene_id'],'ENSG1')
            self.assertEqual(result['scores'][0]['track']['name'],'brain cortex')
            self.assertEqual(result['scores'][0]['quantile_score'],.7)
            self.assertNotIn('FAKE_SECRET',(o/'query_status.tsv').read_text());self.assertNotIn('FAKE_SECRET',(o/'run_config.json').read_text())
            subprocess.run([sys.executable,str(D/'integrate.py'),'--evidence',str(e),'--run',str(o),'--out',str(p/'joined')],check=True,capture_output=True)
            joined=rows(p/'joined'/'variant_molecular_evidence.tsv')
            self.assertEqual([r['atlas_query_status'] for r in joined],['SCORED','QUERY_FAILED','MAPPING_HOLD'])
    def test_jbrowse_six_predictions_and_window_summaries(self):
        combined=OUT/'jbrowse_predictions';receipt=json.loads((combined/'provenance.json').read_text())
        for name,digest in receipt['output_sha256'].items():self.assertEqual(sha(combined/name),digest)
        summaries=rows(combined/'six_region_predictions.tsv');self.assertEqual(len(summaries),6)
        self.assertEqual({s['candidate_locus_id'] for s in summaries},{s['candidate_locus_id'] for s in self.r})
        evid=rows(OUT/'avi_mirror'/'variant_evidence_avi.tsv')
        for record in summaries:
            folder=OUT/record['prediction_folder'];prov=json.loads((folder/'provenance.json').read_text())
            for name,digest in prov['output_sha256'].items():self.assertEqual(sha(folder/name),digest)
            self.assertNotIn('https://',json.dumps(prov['tracks']))
            self.assertNotIn('token',prov)
            variant=prov['request']['variant']
            self.assertEqual(f"{variant['chromosome']}:{variant['position']}:{variant['reference_bases']}>{variant['alternate_bases']}",record['variant_grch38'])
            group=[r for r in evid if r['candidate_locus_id']==record['candidate_locus_id']]
            expected=sorted(group,key=lambda x:(-int(x['avi_raw_e5']),x['atlas_variant']))[0]
            self.assertEqual(expected['atlas_variant'],record['variant_grch38'])
            w=rows(folder/'track_window.tsv');tracks=rows(folder/'track_summary.tsv')
            self.assertEqual(len(w),8004);self.assertEqual(len(tracks),4)
            for t in tracks:
                series=[x for x in w if x['output_type']==t['output_type'] and x['track_ordinal']==t['track_ordinal']]
                self.assertEqual(len(series),2001);self.assertEqual([int(x['offset_from_variant_bp']) for x in series],list(range(-1000,1001)))
                mid=[x for x in series if -250<=int(x['offset_from_variant_bp'])<=250]
                ref=sum(float(x['reference']) for x in mid);alt=sum(float(x['alternate']) for x in mid)
                self.assertAlmostEqual(ref,float(t['ref_sum_501bp']),places=10);self.assertAlmostEqual(alt,float(t['alt_sum_501bp']),places=10)
                self.assertAlmostEqual(math.log2((alt+1)/(ref+1)),float(t['local_log2_ratio_plus_one']),places=10)
if __name__=='__main__':unittest.main(verbosity=2)

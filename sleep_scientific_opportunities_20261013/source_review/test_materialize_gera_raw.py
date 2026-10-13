"""Small SYNTHETIC invariants; no scientific source/reference body or fit."""
import csv
from fractions import Fraction
import gzip
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).with_name('materialize_gera_raw.py')
spec = importlib.util.spec_from_file_location('gera_adapter_fixture', SCRIPT)
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)
SCRATCH = Path('/Volumes/Extreme SSD/sleep_scientific_opportunities_20261013/source_review')


def source_row(bp, **changes):
    row = dict.fromkeys(adapter.HEADER, '')
    row.update(chromosome='1', base_pair_location=str(bp), rs_id_all='provenance_not_an_rsID',
               info_all='0.95', other_allele='C', effect_allele='A', all_maf='0.2',
               all_total='56637', standard_error='0.1', beta='0.2', p_value='0.9')
    row.update(changes)
    return row


class SyntheticScientificInvariants(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='GERA_SYNTHETIC_ONLY_', dir=SCRATCH)
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def run_fixture(self, refs, rows, allele_rows=None, folder='out'):
        source = self.root/'synthetic_source.tsv'
        reference = self.root/'synthetic_reference.tsv.gz'
        alleles = self.root/'synthetic_alleles.tsv'
        with source.open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=adapter.HEADER, delimiter='\t', lineterminator='\n')
            w.writeheader()
            w.writerows(rows)
        with gzip.open(reference, 'wt', newline='') as f:
            w = csv.writer(f, delimiter='\t', lineterminator='\n')
            w.writerow(('SNP','CHR','BP'))
            w.writerows((snp, chromosome, bp) for snp,chromosome,bp,a1,a2 in refs)
        with alleles.open('w', newline='') as f:
            w = csv.writer(f, delimiter='\t', lineterminator='\n')
            w.writerow(('SNP','A1','A2'))
            w.writerows(allele_rows if allele_rows is not None else
                        dict((snp,(snp,a1,a2)) for snp,chr,bp,a1,a2 in refs).values())
        data = source.read_bytes()
        identity = {'bytes': len(data), 'md5': hashlib.md5(data).hexdigest(),
                    'sha256': hashlib.sha256(data).hexdigest()}
        qc = adapter.materialize(source, reference, alleles, self.root/folder, identity, minimum_reference=1)
        output = []
        if qc['sumstats_path']:
            with gzip.open(qc['sumstats_path'], 'rt') as f:
                output = list(csv.DictReader(f, delimiter='\t'))
        self.assertEqual(qc['input_rows'], qc['reconciled_rows'])
        return qc, output

    def test_unique_coordinate_all_orientations_and_total_N_not_Neff(self):
        refs = [('rs1',1,100,'A','G'), ('rs9',1,100,'A','C'),
                ('rs2',1,101,'A','C'), ('rs3',1,102,'A','C'), ('rs4',1,103,'A','C')]
        rows = [source_row(100, rs_id_all='rs1;wrong_first_label'),
                source_row(101,effect_allele='C',other_allele='A'),
                source_row(102,effect_allele='T',other_allele='G'),
                source_row(103,effect_allele='G',other_allele='T',all_total='56637.0')]
        qc, out = self.run_fixture(refs, rows)
        self.assertEqual([x['SNP'] for x in out], ['rs9','rs2','rs3','rs4'])
        self.assertEqual([float(x['Z']) for x in out], [2,-2,2,-2])
        self.assertEqual({(x['A1'],x['A2']) for x in out}, {('A','C')})
        self.assertEqual({int(x['N']) for x in out}, {56637})
        self.assertNotEqual(float(out[0]['N']),14773.761322103925)
        self.assertEqual(qc['reference_counts']['coordinate_collision_groups'],1)
        self.assertEqual(qc['orientation_counts_mapped_input_rows'],
                         {'same':1,'swap':1,'complement':1,'complement_swap':1})
        self.assertGreater(qc['diagnostic_flags']['max_abs_source_P_minus_Wald_P'],0.8)

    def test_multiple_compatible_rsids_hold_no_native_input(self):
        qc,out = self.run_fixture([('rs1',1,100,'A','C'),('rs2',1,100,'A','C')], [source_row(100)])
        self.assertFalse(out)
        self.assertIsNone(qc['sumstats_path'])
        self.assertEqual(qc['status'],'HOLD_NO_NATIVE_INPUT')
        self.assertEqual(qc['exclusive_exclusions']['multiple_compatible_rsids'],1)

    def test_invalid_SE_or_N_counterpart_still_excludes_complete_SNP_group(self):
        refs = [('rs1',1,100,'A','C'),('rs2',1,200,'A','C'),
                ('rs3',1,300,'A','C'),('rs4',1,400,'A','C')]
        first=source_row(100)
        rows=[first,first.copy(),source_row(100,standard_error='0'),source_row(200),
              source_row(200,all_total='56637.000000000000000000001'),
              source_row(300),source_row(400,all_total='56638')]
        qc,out=self.run_fixture(refs,rows)
        self.assertEqual([x['SNP'] for x in out],['rs3'])
        self.assertEqual(qc['duplicate_copies'],1)
        self.assertEqual(qc['conflicting_groups'],2)
        self.assertEqual(qc['conflicting_unique_records'],4)
        self.assertEqual(qc['conflicting_source_rows_diagnostic'],5)
        self.assertEqual(qc['exclusive_exclusions']['N_above_56637'],1)
        self.assertEqual(qc['diagnostic_flags']['N_fractional'],1)
        self.assertEqual(qc['diagnostic_flags']['SE_not_positive'],1)
        self.assertEqual(qc['N_quantile_population_unique_nonconflicting_eligible_mapped'],1)

    def test_quantile_copy_weights_excluded_and_equality_retained(self):
        refs=[('rs'+str(i+1),1,100+i,'A','C') for i in range(5)]
        low=source_row(100,all_total='99')
        rows=[low.copy() for _ in range(21)]+[source_row(101,all_total='1e2')]+[
            source_row(100+i,all_total='150') for i in (2,3,4)]
        qc,out=self.run_fixture(refs,rows)
        self.assertEqual(qc['R7_exact_m_j_r_Q'],{'m':5,'j':3,'r':6,'Q':1500})
        self.assertEqual(qc['duplicate_copies'],20)
        self.assertEqual(qc['N_filtered_unique_records'],1)
        self.assertEqual(qc['N_equal_threshold'],1)
        self.assertEqual(qc['N_above_threshold'],3)
        self.assertEqual([x['SNP'] for x in out],['rs2','rs3','rs4','rs5'])
        self.assertEqual(adapter.r7_histogram({10:1,20:1,30:1,40:1,50:1},Fraction(9,10)),46)

    def test_reference_duplicate_keys_collapse_identical_and_taint_all_conflicting_coordinates(self):
        refs=[('rs1',1,100,'A','C'),('rs1',1,200,'A','C'),
              ('rs3',1,300,'A','C'),('rs3',1,300,'A','C')]
        qc,out=self.run_fixture(refs,[source_row(100),source_row(200),source_row(300)],
                                allele_rows=[('rs1','A','C'),('rs3','A','C'),('rs3','A','C')],
                                folder='ref_duplicate')
        self.assertEqual([x['SNP'] for x in out],['rs3'])
        self.assertEqual(qc['exclusive_exclusions']['tainted_reference_coordinate'],2)
        self.assertEqual(qc['reference_counts']['identical_reference_definition_copies'],1)
        self.assertEqual(qc['reference_counts']['identical_allele_definition_copies'],1)
        self.assertEqual(qc['reference_counts']['conflicting_reference_SNP_groups'],1)
        # An unrelated good SNP survives an allele-list conflict at another SNP.
        refs=[('rs1',1,100,'A','C'),('rs4',1,400,'A','C')]
        qc,out=self.run_fixture(refs,[source_row(100),source_row(400)],
                               allele_rows=[('rs1','A','C'),('rs1','A','G'),('rs4','A','C')],
                               folder='allele_duplicate')
        self.assertEqual([x['SNP'] for x in out],['rs4'])
        self.assertEqual(qc['reference_counts']['conflicting_allele_SNP_groups'],1)
        self.assertEqual(qc['exclusive_exclusions']['tainted_reference_coordinate'],1)

    def test_strict_filters_MHC_boundaries_and_unclipped_INFO_anomaly(self):
        refs=[('rs1',1,100,'A','C'),('rs2',1,101,'A','C'),('rs3',1,102,'A','C'),
              ('rs4',6,25000000,'A','C'),('rs5',6,34000000,'A','C'),
              ('rs6',6,24999999,'A','C'),('rs7',6,34000001,'A','C')]
        rows=[source_row(100,info_all='0.9'),source_row(101,all_maf='0.01'),
              source_row(102,info_all='1.2',all_maf='0.5')]+[
                  source_row(bp,chromosome='6') for bp in (25000000,34000000,24999999,34000001)]
        qc,out=self.run_fixture(refs,rows)
        self.assertEqual([x['SNP'] for x in out],['rs3','rs6','rs7'])
        self.assertEqual(qc['exclusive_exclusions']['extended_MHC37'],2)
        self.assertEqual(qc['exclusive_exclusions']['INFO_not_gt_0.9'],1)
        self.assertEqual(qc['exclusive_exclusions']['MAF_not_in_0.01_0.5'],1)
        self.assertEqual(qc['diagnostic_flags']['INFO_outside_0_1'],1)


if __name__ == '__main__':
    unittest.main(verbosity=2)

"""Structural dry run for the final join while real secondary analysis is pending."""
from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / 'build_brain6_final_package_v1.py'


class FinalPackageBuilderTest(unittest.TestCase):
    def test_complete_join_with_synthetic_secondary_statuses(self) -> None:
        spec = importlib.util.spec_from_file_location('brain6_final_builder_test', SCRIPT)
        assert spec and spec.loader
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory(prefix='brain6_final_join_') as directory:
            root = Path(directory)
            alternative = root / 'alternative'
            alternative.mkdir()
            base = module.table(module.BOUNDED / 'candidate_evidence_25.tsv')
            regions = module.table(module.BOUNDED / 'region_evidence_20.tsv')
            with (alternative / 'candidate_crosswalk.tsv').open('w', newline='') as stream:
                writer = csv.DictWriter(stream, delimiter='\t',
                                        fieldnames=['candidate_locus_id', 'pair_id', 'status',
                                                   'source_eligibility', 'best_p', 'covariance_sign'],
                                        lineterminator='\n')
                writer.writeheader()
                for row in base:
                    writer.writerow({'candidate_locus_id': row['candidate_locus_id'],
                                     'pair_id': row['pair_id'],
                                     'status': ('METHOD_INAPPLICABLE' if row['pair_id'].startswith('longsleep__')
                                                else 'NOT_TESTED_INPUT_QC'),
                                     'source_eligibility': ('METHOD_INAPPLICABLE' if row['pair_id'].startswith('longsleep__')
                                                            else 'ELIGIBLE_SECONDARY'), 'best_p': 'NA',
                                     'covariance_sign': 'NA'})
            with (alternative / 'region_crosswalk.tsv').open('w', newline='') as stream:
                writer = csv.DictWriter(stream, delimiter='\t', fieldnames=['region_group'], lineterminator='\n')
                writer.writeheader()
                for row in regions:
                    writer.writerow({'region_group': row['geographic_region_grch37']})
            pair_status = {pair: ('METHOD_INAPPLICABLE' if pair.startswith('longsleep__') else 'ELIGIBLE_SECONDARY')
                           for pair in sorted({row['pair_id'] for row in base})}
            (alternative / 'source_eligibility_lock_v1.json').write_text(json.dumps({
                'pair_status': pair_status, 'family_denominator_pairs': 5,
                'blocks_per_pair': 1693, 'alpha_block': 0.05 / 8465}))
            (alternative / 'run_provenance.json').write_text('{}')
            with (alternative / 'genome_wide_blocks.tsv').open('w', newline='') as stream:
                writer = csv.DictWriter(stream, delimiter='\t',
                                        fieldnames=['pair_id', 'analysis_status'], lineterminator='\n')
                writer.writeheader()
                for pair in pair_status:
                    for _ in range(1693):
                        writer.writerow({'pair_id': pair,
                                         'analysis_status': ('METHOD_INAPPLICABLE' if pair.startswith('longsleep__')
                                                             else 'NOT_TESTED_LOW_SNP')})
            module.ALT = alternative
            module.OUT = root / 'package'
            module.main()
            self.assertEqual(len(module.table(module.OUT / 'BRAIN6_FINAL_CANDIDATES.tsv')), 25)
            self.assertEqual(len(module.table(module.OUT / 'BRAIN6_FINAL_REGIONS.tsv')), 20)
            self.assertEqual(len(module.table(module.OUT / 'BRAIN6_FINAL_LAVA.tsv')), 17465)
            self.assertEqual(len(module.table(module.OUT / 'BRAIN6_FINAL_GLOBAL.tsv')), 72)
            self.assertEqual(len(module.table(module.OUT / 'BRAIN6_FINAL_ALT_LOCAL_BLOCKS.tsv')), 8465)
            self.assertEqual(len(module.table(module.OUT / 'BRAIN6_FINAL_TISSUE_CELLTYPE.tsv')), 79)
            self.assertEqual(len(module.table(module.OUT / 'BRAIN6_FINAL_PATHWAYS.tsv')), 1700)
            self.assertFalse(any(row['final_evidence_class'] == 'EXPLORATORY_MULTIOMIC_SUPPORTED'
                                 for row in module.table(module.OUT / 'BRAIN6_FINAL_CANDIDATES.tsv')))
            audit_script = SCRIPT.with_name('audit_brain6_final_package_v1.py')
            audit_spec = importlib.util.spec_from_file_location('brain6_final_audit_test', audit_script)
            assert audit_spec and audit_spec.loader
            audit_module = importlib.util.module_from_spec(audit_spec)
            audit_spec.loader.exec_module(audit_module)
            audit_module.PACKAGE = module.OUT
            for name in audit_module.DOCUMENTS:
                content = '# Synthetic document for structural audit\n'
                if name == 'BRAIN6_FINAL_RESULTS.md':
                    content += ('17,465 13,745 3,720 873 25 pair-specific 20 overlapping '
                                '47 `TESTED` 41 `NOT_RUN` 3,304 5,079 24 model-conditional '
                                'maximum PP.H4 was **0.0824** 0/25 54 GTEx 1,680 Reactome '
                                'A `PRIMARY_LAVA_CONFIRMED`: 0 E `UNSUPPORTED`: 25\n')
                if name == 'BRAIN6_FINAL_CLAIMS.md':
                    content += ('PRIMARY_CONFIRMATORY_SUPPORTED SECONDARY_VALIDATION_SUPPORTED '
                                'EXPLORATORY_SUPPORTED DESCRIPTIVE_ONLY UNSUPPORTED_DO_NOT_CLAIM\n')
                (module.OUT / name).write_text(content)
            self.assertEqual(audit_module.audit()['status'], 'PASS')


if __name__ == '__main__':
    unittest.main()

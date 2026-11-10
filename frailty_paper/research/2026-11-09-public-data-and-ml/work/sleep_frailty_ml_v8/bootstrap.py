"""Acquire pinned small-model sources and fresh physical designs, no outcomes."""
from v8_common import *

COMMIT = 'a20dbc7d0b08daae3af1aebcfb8a71d9f198ab5d'


def main():
    assert_previous_intact()
    plan_path = OUT/'initial_experiment_plan.json'
    if not plan_path.exists():
        write_json(plan_path, {
            'purpose': 'Improve measured magnitude, nearby-variant and cell-specific transfer; no person-level disease prediction.',
            'history': 'V7 outside failures are revealed benchmarks. They cannot be treated as untouched tests for V8.',
            'base_model': 'Author-distributed MPRA-LegNet WTC11 example checkpoint, pinned repository revision. Exact checkpoint fold must be established or reported unknown.',
            'training_candidates': ['Source-verified signed laboratory effects and absolute magnitudes from the existing V7 development sources', 'Frozen MPRA-trained features versus actual partial weight fine-tuning', 'Signed allele-difference prediction plus magnitude-focused supervised prediction', 'Matched measured cell differences and nearby magnitude differences'],
            'development_only_selection': 'Keep linked regions, identical DNA and all repeated cell outcomes together. Evaluate entire development studies left out in turn. Lock the recipe, training rows, confidence rules and all predictions before opening fresh-study effect performance.',
            'fresh_study_candidate': 'Matoba et al., GSE289858. Physical DNA, coordinate/strand, original effect conventions, duplicates, donor context and training/pretraining overlap must be audited first. Availability does not establish eligibility.',
            'AlphaGenome_input': 'Only designated permissible composite AVI; matched with/without experiment. No detailed API-output teacher labels. Retain only if added value survives a fresh-study comparison and no material local/cell ranking harm.',
            'success_requirements': ['Meaningful lower measured-magnitude error than training-cell medians on whole-study transfer', 'Variant-specific cell differences improve over both median-shift and zero-difference comparators', 'Nearby ordering improves on measured comparisons', 'Coverage, interval width and nonzero answer availability reported jointly; no all-abstention success claim', 'Independent fresh-study test after recipe and predictions are locked'],
            'no_success_claim_from': ['Checkpoint loading', 'Nonzero gradients', 'A completed training job', 'Training-set accuracy', 'Copying AlphaGenome outputs', 'Withholding all predictions'],
            'original_evidence_preserved': True,
        })
    receipts = []
    for name in ['model.py', 'dataset.py', 'utils.py', 'LICENSE', 'models/example.cfg', 'models/example.cpkt', 'datasets/original/WTC11.tsv']:
        target = SOURCES/('legnet_'+name.replace('/', '_'))
        receipts.append(fetch(f'https://raw.githubusercontent.com/autosome-ru/human_legnet/{COMMIT}/{name}', target, 30_000_000))
        print('SOURCE ACQUIRED', name, target.stat().st_size, flush=True)
    url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE289nnn/GSE289858/suppl/GSE289858_variant_barcode_statistics.txt.gz'
    receipts.append(fetch(url, SOURCES/'GSE289858_variant_barcode_statistics.txt.gz', 25_000_000, b'\x1f\x8b'))
    receipts.append(fetch('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE289nnn/GSE289858/soft/GSE289858_family.soft.gz', SOURCES/'GSE289858_family.soft.gz', 5_000_000, b'\x1f\x8b'))
    write_json(OUT/'source_acquisition_receipt.json', {'legnet_repository_commit': COMMIT, 'sources': receipts, 'Matoba_model_performance_calculated': False, 'source_files_modified': False})
    print('V8 SOURCES READY; FRESH OUTCOME PERFORMANCE STILL UNTOUCHED', flush=True)


if __name__ == '__main__':
    main()

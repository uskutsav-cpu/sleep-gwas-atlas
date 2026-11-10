"""Correct the descriptive label, preserving all original coefficients/evidence."""
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/sleep_frailty_ml_neuronal_activity'
OLD = 'Original construct-level logFC from BCalm topTreat versus scrambled controls; not an allelic contrast or phenotype risk.'
NEW = ('Original unshifted BCalm fitted construct log2 RNA/DNA activity coefficient. '
       'The scrambled-control choice determines TREAT statistics; the released logFC is not the difference from the scrambled-control mean. '
       'Not an allelic contrast or phenotype risk.')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    path = OUT / 'additional_neuronal_activity_inventory_receipt.json'
    before = path.read_bytes()
    data = json.loads(before)
    archive = OUT / 'inventory_receipt_before_activity_target_description_correction.json'
    if data['source_target'] == NEW:
        assert archive.exists()
        return
    assert data['source_target'] == OLD
    paper = OUT / 'BCalm_original_2025_methods.html'
    assert b'original unshifted logratios' in paper.read_bytes()
    code_paths = [OUT / 'BCalm_v0_99_analyze.R', OUT / 'BCalm_current_analyze.R']
    for p in code_paths:
        code = p.read_text()
        assert 'coefficients <- as.matrix(fit$coefficients - neg_mean)' in code
        assert not re.search(r'fit\$coefficients\s*(<-|=)', code)
    cli = OUT / 'perform_bcalm_elements_cli.R'
    assert 'topTreat(testing_element_scrambled, coef = 1, number = Inf)' in cli.read_text()
    if archive.exists():
        assert archive.read_bytes() == before
    else:
        archive.write_bytes(before)
    correction_path = OUT / 'activity_target_interpretation_correction.json'
    correction = {'corrected': True, 'old_ambiguous_description': OLD, 'correct_description': NEW,
        'why': 'BCalm fits element intercepts on the log2 RNA/DNA scale. mpra_treat centers a local coefficient copy to calculate control-based statistics, while returning unshifted coefficients. The source CLI writes topTreat logFC values.',
        'primary_methods': {'url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC11827149/',
                            'local_HTML_sha256': sha(paper), 'section': 'Materials and methods, Algorithm'},
        'historical_source': {'url': 'https://github.com/kircherlab/BCalm/blob/cf94eedd5d3a96a2e8794b75aeaeeaaf034a47a8/R/analyze.R',
                              'tag': 'v0.99.0', 'sha256': sha(code_paths[0]),
                              'exact_software_checkout_used_by_neuron_authors_established': False},
        'current_source': {'url': 'https://github.com/kircherlab/BCalm/blob/55fecd0e99c1a2eb633165dd90be23277138aa95/R/analyze.R',
                           'sha256': sha(code_paths[1])},
        'original_neuron_CLI_sha256': sha(cli), 'archived_original_receipt_sha256': sha(archive),
        'source_metadata_or_original_coefficients_changed': False, 'training_arrays_or_running_comparison_changed': False,
        'statistical_significance_or_t_statistics_used_to_select_new_training_examples': False,
        'raw_counts_or_historical_BCalm_runtime_refitted': False, 'code_sha256': sha(Path(__file__))}
    correction_path.write_text(json.dumps(correction, indent=2) + '\n')
    data['source_target'] = NEW
    data['source_target_description_correction'] = {'receipt': str(correction_path.relative_to(ROOT)),
                                                  'sha256': sha(correction_path),
                                                  'original_receipt_archive': str(archive.relative_to(ROOT))}
    path.write_text(json.dumps(data, indent=2) + '\n')
    print(json.dumps(correction), flush=True)


if __name__ == '__main__':
    main()

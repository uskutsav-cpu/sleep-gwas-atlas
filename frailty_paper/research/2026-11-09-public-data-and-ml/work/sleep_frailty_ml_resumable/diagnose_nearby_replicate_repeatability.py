"""Nearby-variant ordering across real assay replicates, not a model ceiling.

Recover original Kosicki RNA/DNA readouts, verify their full current-development
identity and compare disjoint replicate sets. No calibration-region outcome,
reserved-study outcome, model prediction, significance filter or fit is used.
"""
import csv
import datetime
import gzip
import itertools
import json
from pathlib import Path

import numpy as np

from readonly_native_inputs import ROOT, sha
from audit_saved_pair_metrics import correlation

OUT = ROOT / 'outputs/sleep_frailty_ml_resumable'
BASE = Path('/Volumes/Extreme SSD/codex-sleep-frailty-work')
V6, V7, V8 = [BASE / ('sleep_frailty_ml_v' + str(v)) for v in [6, 7, 8]]


def compare(a, b):
    keep = (a != 0) & (b != 0)
    result = {'pairs': len(a), 'pairs_without_an_exact_zero_contrast': int(keep.sum()),
              'excluded_exact_ties': int((~keep).sum())}
    if not keep.any():
        return result
    aa, bb = a[keep], b[keep]
    pa, pb = float((aa > 0).mean()), float((bb > 0).mean())
    agreement = float((np.sign(aa) == np.sign(bb)).mean())
    chance = pa * pb + (1 - pa) * (1 - pb)
    result.update(ordering_agreement=agreement,
                  positive_contrast_fractions=[pa, pb],
                  marginal_sign_expected_agreement=chance,
                  chance_adjusted_sign_agreement=(agreement - chance) / (1 - chance) if chance < 1 else None,
                  contrast_Spearman=correlation(aa, bb),
                  contrast_MAE=float(abs(aa - bb).mean()),
                  zero_contrast_MAE_for_second_set=float(abs(bb).mean()))
    return result


def main():
    destination = OUT / 'nearby_assay_replicate_repeatability.json'
    assert not destination.exists()
    data = ROOT / 'outputs/sleep_frailty_ml_next/coverage_training_arrays.npz'
    pairs_path = data.with_name('coverage_measured_pairs.npz')
    plan = json.loads((OUT / 'finite_stable_native_plan.json').read_text())
    expansion = json.loads((data.parent / 'coverage_expansion_receipt.json').read_text())
    assert sha(data) == plan['data_sha256'] and sha(pairs_path) == expansion['pairs_sha256']
    keys = ['bank_index', 'cell', 'line', 'study', 'experiment', 'group', 'split',
            'effect', 'signed_effect', 'original_V8_observation_index', 'study_names']
    with np.load(data, allow_pickle=False) as archive:
        d = {key: archive[key] for key in keys}
    sid = d['study_names'].tolist().index('Kosicki2025')
    rows = np.flatnonzero((d['study'] == sid) & (d['split'] != 2))
    assert len(np.unique(d['bank_index'][rows])) == len(rows)
    bank_to_row = dict(zip(map(int, d['bank_index'][rows]), map(int, rows)))
    metadata = json.loads((OUT / 'development/frozen_reader_Kosicki2025_42_started.json').read_text())
    sequence_path = V7 / 'sequence_bank.csv.gz'
    source_hashes = {entry['path']: entry['sha256'] for entry in metadata['physical_source_receipt']['sources']}
    assert sha(sequence_path) == source_hashes[str(sequence_path)]
    design_to_row, positions, chromosomes = {}, {}, {}
    with gzip.open(sequence_path, 'rt') as handle:
        for record in csv.DictReader(handle):
            bank = int(record['bank_id'])
            if bank in bank_to_row:
                assert record['context_feature_source'] == 'kosicki'
                design = int(record['context_feature_index'])
                assert design not in design_to_row
                design_to_row[design] = bank_to_row[bank]
                positions[bank], chromosomes[bank] = int(record['position_hg38']), record['chromosome']
    assert set(design_to_row.values()) == set(rows)
    names_by_row = {}
    reconstructed = V6 / 'outputs/kosicki_reconstructed_continuous_effects.csv'
    prior = json.loads((V6 / 'outputs/kosicki_continuous_recovery_audit.json').read_text())
    assert sha(reconstructed) == prior['reconstructed_effects_sha256']
    with reconstructed.open() as handle:
        for record in csv.DictReader(handle):
            row = design_to_row.get(int(record['design_index']))
            if row is not None:
                assert row not in names_by_row
                names_by_row[row] = (record['reference_name'], record['name'])
    assert set(names_by_row) == set(rows)
    names = {name for pair in names_by_row.values() for name in pair}
    original = V6 / 'reserved_outcome_source/published_replicate_readouts.csv'
    assert sha(original) == prior['replicate_export_sha256']
    measurements = {}
    with original.open() as handle:
        for record in csv.DictReader(handle):
            name = record['name']
            if name in names:
                key = (name, int(record['rep']))
                assert key not in measurements
                measurements[key] = float(record['log2.ratio'])
    differences = np.array([[measurements[(names_by_row[int(row)][1], rep)] -
                             measurements[(names_by_row[int(row)][0], rep)] for rep in [1, 2, 3]]
                            for row in rows], dtype='float64')
    assert np.isfinite(differences).all()
    magnitude_error = float(abs(abs(differences.mean(1)) - d['effect'][rows]).max())
    signed_error = float(abs(differences.mean(1) - d['signed_effect'][rows]).max())
    assert max(magnitude_error, signed_error) < 2e-6
    cached = V8 / 'development_Kosicki_actual_log2_replicate_differences.npz'
    with np.load(cached, allow_pickle=False) as archive:
        old = archive['observation_index']; values = archive['log2_difference']
        assert len(np.unique(old)) == len(old)
        mapping = {int(index): i for i, index in enumerate(old)}
        selected = np.array([mapping[int(index)] for index in d['original_V8_observation_index'][rows]])
        cache_error = float(abs(values[selected] - differences).max())
    assert cache_error == 0
    inverse = np.full(len(d['effect']), -1, dtype='int32'); inverse[rows] = np.arange(len(rows))
    with np.load(pairs_path, allow_pickle=False) as archive:
        pairs = archive['local']
    pairs = pairs[np.all(inverse[pairs] >= 0, axis=1)]
    a, b = pairs.T
    for key in ['cell', 'line', 'study', 'experiment', 'group', 'split']:
        assert np.all(d[key][a] == d[key][b])
    assert not np.any(d['split'][pairs] == 2)
    for ia, ib in zip(d['bank_index'][a], d['bank_index'][b]):
        assert ia != ib and chromosomes[int(ia)] == chromosomes[int(ib)]
        assert abs(positions[int(ia)] - positions[int(ib)]) <= 25
    one = abs(differences)
    pair_delta = one[inverse[a]] - one[inverse[b]]
    pairwise = {f'replicate_{i+1}_versus_{j+1}': compare(pair_delta[:, i], pair_delta[:, j])
                for i, j in itertools.combinations(range(3), 2)}
    disjoint = {}
    for i in range(3):
        other = [j for j in range(3) if j != i]
        magnitude = abs(differences[:, other].mean(1))  # Mean signed difference, then absolute.
        contrast = magnitude[inverse[a]] - magnitude[inverse[b]]
        disjoint[f'replicate_{i+1}_versus_mean_of_others'] = compare(pair_delta[:, i], contrast)
    result = {'completed': True, 'time_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'source': 'Kosicki2025', 'actual_noncalibration_DNA_observations': len(rows),
              'actual_local_pairs': len(pairs), 'connected_genomic_pair_groups': len(np.unique(d['group'][a])),
              'all_three_original_assay_replicate_vectors_reconstructed': True,
              'maximum_current_magnitude_identity_error': magnitude_error,
              'maximum_current_signed_identity_error': signed_error,
              'maximum_prior_replicate_cache_identity_error': cache_error,
              'pairwise_single_replicate_ordering': pairwise,
              'disjoint_one_versus_two_replicate_ordering': disjoint,
              'no_effect_size_significance_or_SE_filter': True,
              'measurement': 'Per replicate ALT-minus-REF log2 RNA/DNA activity. Magnitude uses abs(mean signed differences), not mean absolute differences.',
              'scope': 'Experimental repeatability diagnostic on the established development corpus. Shared controls, loci and assay preparation mean pairs are not independent donor replicates. This is not a formal prediction ceiling, new-study validation, uncertainty calibration or proof that measurement noise caused model failure.',
              'current_models_targets_partitions_and_training_unchanged': True,
              'model_predictions_loaded_or_selected': False, 'reserved_new_study_outcomes_read': False,
              'original_calibration_outcomes_used': False, 'scientific_goal_complete': False,
              'source_sha256': {str(path): sha(path) for path in [data, pairs_path, sequence_path, reconstructed, original, cached]},
              'code_sha256': {path.name: sha(path) for path in [Path(__file__), Path(__file__).with_name('audit_saved_pair_metrics.py')]}}
    destination.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in ['actual_noncalibration_DNA_observations', 'actual_local_pairs', 'connected_genomic_pair_groups', 'maximum_prior_replicate_cache_identity_error', 'pairwise_single_replicate_ordering', 'disjoint_one_versus_two_replicate_ordering']}), flush=True)


if __name__ == '__main__':
    main()

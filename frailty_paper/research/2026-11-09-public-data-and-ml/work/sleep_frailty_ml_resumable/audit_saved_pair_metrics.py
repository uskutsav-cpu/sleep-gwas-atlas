"""Independent measured-pair replay, including exact-line comparisons."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/sleep_frailty_ml_resumable'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def correlation(x, y):
    def ranks(values):
        order = np.argsort(values, kind='stable')
        sorted_values = values[order]
        edges = np.r_[0, np.flatnonzero(sorted_values[1:] != sorted_values[:-1]) + 1, len(values)]
        means = (edges[:-1] + edges[1:] - 1) / 2.
        result = np.empty(len(values), dtype='float64')
        result[order] = np.repeat(means, np.diff(edges))
        return result - result.mean()
    a, b = ranks(x), ranks(y)
    denominator = np.sqrt(np.dot(a, a) * np.dot(b, b))
    return float(np.dot(a, b) / denominator) if denominator else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--recipe', default='frozen_reader')
    parser.add_argument('--source', default='none')
    parser.add_argument('--epoch', type=int, default=3)
    args = parser.parse_args()
    stem = f'{args.recipe}_{args.source}_42'
    folder = OUT / 'development'
    metadata = json.loads((folder / (stem + '_started.json')).read_text())
    receipt = json.loads((folder / f'{stem}_epoch{args.epoch}_receipt.json').read_text())
    path = folder / f'{stem}_epoch{args.epoch}_predictions.npz'
    assert sha(path) == receipt['prediction_sha256']
    plan = json.loads((OUT / 'finite_stable_native_plan.json').read_text())
    data = ROOT / 'outputs/sleep_frailty_ml_next/coverage_training_arrays.npz'
    pairs_path = ROOT / 'outputs/sleep_frailty_ml_next/coverage_measured_pairs.npz'
    assert sha(data) == plan['data_sha256']
    data_receipt = json.loads((ROOT / 'outputs/sleep_frailty_ml_next/coverage_expansion_receipt.json').read_text())
    assert sha(pairs_path) == data_receipt['pairs_sha256']
    with np.load(data, allow_pickle=False) as z:
        d = {k: z[k] for k in ['effect', 'cell', 'line', 'study', 'group', 'split', 'bank_index', 'study_names', 'cell_names', 'line_names'] if k in z.files}
    with np.load(path, allow_pickle=False) as z:
        ids, prediction = z['observation_index'], z['magnitude']
    assert hashlib.sha256(ids.tobytes()).hexdigest() == metadata['prediction_index_sha256']
    train = d['split'] == 0
    if args.source != 'none':
        sid = d['study_names'].tolist().index(args.source)
        excluded = np.unique(d['group'][d['study'] == sid])
        train &= (d['study'] != sid) & ~np.isin(d['group'], excluded)
    assert hashlib.sha256(np.flatnonzero(train).tobytes()).hexdigest() == metadata['training_index_sha256']
    cell_medians = np.full(len(metadata['cell_medians']), np.median(d['effect'][train]), dtype='float32')
    for c in np.unique(d['cell'][train]):
        cell_medians[c] = np.median(d['effect'][train & (d['cell'] == c)])
    assert np.array_equal(cell_medians, np.array(metadata['cell_medians'], dtype='float32'))
    line_medians = cell_medians[np.asarray(metadata['line_parent'])].copy()
    for line in np.unique(d['line'][train]):
        line_medians[line] = np.median(d['effect'][train & (d['line'] == line)])
    assert np.array_equal(line_medians, np.array(metadata['line_medians'], dtype='float32'))
    n = int(d['bank_index'].max()) + 1
    positions, chromosomes = np.full(n, -1, dtype='int64'), np.full(n, -1, dtype='int16')
    chrom_ids = {'chr' + str(i): i for i in range(1, 23)} | {'chrX': 23, 'chrY': 24}
    for entry in metadata['physical_source_receipt']['sources']:
        assert sha(Path(entry['path'])) == entry['sha256']
        with gzip.open(entry['path'], 'rt') as f:
            for row in csv.DictReader(f):
                i = int(row['bank_id'])
                if i < n:
                    positions[i], chromosomes[i] = int(row['position_hg38']), chrom_ids[row['chromosome']]
    inverse = np.full(len(d['effect']), -1, dtype='int32')
    inverse[ids] = np.arange(len(ids))
    curve = json.loads((folder / (stem + '_curve.json')).read_text())[args.epoch - 1]
    results, errors = {}, []
    with np.load(pairs_path, allow_pickle=False) as z:
        for label in z.files:
            pair = z[label]
            pair = pair[np.all(inverse[pair] >= 0, axis=1)]
            if not len(pair):
                results[label] = {'pairs': 0}
                assert curve['comparisons'][label]['pairs'] == 0
                continue
            a, b = pair[:, 0], pair[:, 1]
            assert np.all(d['split'][a] == d['split'][b])
            assert np.all(d['study'][a] == d['study'][b])
            assert np.all(d['group'][a] == d['group'][b])
            if label == 'local':
                assert np.all(d['line'][a] == d['line'][b])
                ba, bb = d['bank_index'][a], d['bank_index'][b]
                assert np.all(ba != bb)
                assert np.all(chromosomes[ba] == chromosomes[bb]) and np.all(chromosomes[ba] > 0)
                assert np.all(abs(positions[ba] - positions[bb]) <= 25)
            else:
                assert label == 'cross'
                assert np.all(d['bank_index'][a] == d['bank_index'][b])
                assert np.all(d['line'][a] != d['line'][b])
            dy = d['effect'][a] - d['effect'][b]
            dp = prediction[inverse[a]] - prediction[inverse[b]]
            error = float(np.mean(abs(dy - dp), dtype='float64'))
            zero = float(np.mean(abs(dy), dtype='float64'))
            nonzero = dy != 0
            summary = {'pairs': len(pair), 'difference_MAE': error, 'zero_difference_MAE': zero,
                'skill_vs_zero_difference': 1 - error / zero,
                'difference_Spearman': correlation(dy, dp),
                'ordering_accuracy_nonzero_measured_differences': float(np.count_nonzero((np.sign(dy) == np.sign(dp)) & nonzero) / nonzero.sum())}
            for key, value in summary.items():
                errors.append(abs(value - curve['comparisons'][label][key]))
            if label == 'cross':
                mean_cell = cell_medians[d['cell'][a]] - cell_medians[d['cell'][b]]
                cell_error = float(np.mean(abs(dy - mean_cell), dtype='float64'))
                line_diff = line_medians[d['line'][a]] - line_medians[d['line'][b]]
                line_error = float(np.mean(abs(dy - line_diff), dtype='float64'))
                strata = list(zip(d['study'][a], d['cell'][a], d['cell'][b]))
                grouped = {}
                for i, key in enumerate(strata):
                    grouped.setdefault(key, []).append(i)
                values = []
                line_pairs = set()
                for indexes in grouped.values():
                    indexes = np.array(indexes)
                    line_pairs.add(len(set(zip(d['line'][a[indexes]], d['line'][b[indexes]]))))
                    if len(indexes) >= 30 and np.ptp(dy[indexes]) > 0 and np.ptp(dp[indexes]) > 0:
                        values.append(correlation(dy[indexes], dp[indexes]))
                mean_rho = float(np.mean(values))
                errors += [abs(cell_error - curve['comparisons'][label]['training_cell_median_difference_MAE']),
                           abs(mean_rho - curve['comparisons'][label]['macro_within_fixed_study_cell_pair_Spearman'])]
                summary.update(training_cell_median_difference_MAE=cell_error,
                    exact_line_median_difference_MAE=line_error, skill_vs_exact_line_median_difference=1 - error / line_error,
                    macro_within_fixed_study_cell_pair_Spearman=mean_rho, fixed_cell_pair_strata=len(values),
                    exact_line_pair_counts_inside_reported_cell_pair_strata=sorted(line_pairs))
            results[label] = summary
    assert max(errors) < 1e-6, max(errors)
    result = {'passed': True, 'checkpoint': stem + f'_epoch{args.epoch}', 'prediction_sha256': sha(path),
        'actual_pair_identity_and_region_rules_verified': True, 'training_cell_and_line_medians_independently_reconstructed': True,
        'independent_average_tie_rank_correlation_used': True, 'maximum_metric_replay_error': max(errors),
        'comparisons': results, 'new_external_study_or_confidence_success_established': False,
        'scope': 'Saved development predictions and all included measured local/cross-cell pairs. No assertion of donor independence, LD blocks or external biological validation.',
        'code_sha256': sha(Path(__file__))}
    target = OUT / f'{stem}_epoch{args.epoch}_independent_pair_audit.json'
    target.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()

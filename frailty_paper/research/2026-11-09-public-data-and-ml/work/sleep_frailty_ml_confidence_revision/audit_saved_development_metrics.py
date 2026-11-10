"""Independent arithmetic/identity replay of completed confidence experiments."""
import collections
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/sleep_frailty_ml_confidence_revision'
OLD = ROOT / 'outputs/sleep_frailty_ml_next/validation'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    plan_path = OUT / 'selected_confidence_development_plan.json'
    plan = json.loads(plan_path.read_text())
    for name, expected in plan['code_sha256'].items():
        assert sha(Path(__file__).with_name(name)) == expected
    settings_path = OLD / 'final_settings_lock.json'
    assert sha(settings_path) == plan['base_model_settings_sha256']
    settings = json.loads(settings_path.read_text())
    data = ROOT / 'outputs/sleep_frailty_ml_next/coverage_training_arrays.npz'
    assert sha(data) == settings['source_data_sha256']
    with np.load(data, allow_pickle=False) as z:
        y, group, cells, lines = [z[k] for k in ['effect', 'group', 'cell', 'line']]
    audits = {}
    for source in plan['outer_studies']:
        result_path = OUT / f'selected_confidence_{source}_development.json'
        decision_path = OUT / f'selected_confidence_{source}_decisions.npz'
        if not result_path.exists() or not decision_path.exists():
            continue
        result = json.loads(result_path.read_text())
        assert result['plan_sha256'] == sha(plan_path)
        predictions = []
        ids = None
        records = []
        for seed in settings['seeds']:
            record = json.loads((OLD / 'fits' / f'fixed_{source}_sequence_{seed}_record.json').read_text())
            assert record['settings_lock_sha256'] == sha(settings_path)
            assert sha(record['predictions']) == record['predictions_sha256']
            with np.load(record['predictions'], allow_pickle=False) as z:
                current = z['observation_index']
                if ids is None:
                    ids = current
                else:
                    assert np.array_equal(ids, current)
                predictions.append(z['magnitude'])
            records.append(record)
        prediction = np.mean(predictions, axis=0)
        meta = records[0]
        baseline = np.where(np.asarray(meta['supported_lines'])[lines[ids]],
                            np.asarray(meta['line_medians'])[lines[ids]], np.asarray(meta['cell_medians'])[cells[ids]])
        with np.load(decision_path, allow_pickle=False) as z:
            assert np.array_equal(z['observation_index'], ids)
            accepted = z['candidate_accepted']
            width = z['candidate_interval_halfwidth']
            assert np.all(~accepted | z['candidate_eligible'])
            assert np.isfinite(width[accepted]).all() and np.isnan(width[~accepted]).all()
        expected = result['selected_calibration_candidate']
        assert len(ids) == expected['rows'] and int(accepted.sum()) == expected['accepted_rows']
        assert len(set(map(int, group[ids][accepted]))) == expected['accepted_regions']
        checks = {}
        if accepted.any():
            error = abs(y[ids][accepted] - prediction[accepted])
            berror = abs(y[ids][accepted].astype('float64') - baseline[accepted])
            coverage = error <= width[accepted]
            counts = collections.Counter(map(int, group[ids][accepted]))
            weights = np.array([1. / counts[int(g)] for g in group[ids][accepted]])
            row_coverage = sum(map(bool, coverage)) / len(coverage)
            region_coverage = sum(float(w) for w, ok in zip(weights, coverage) if ok) / sum(map(float, weights))
            checks = {'raw_observed_MAE': float(np.mean(error, dtype='float64')),
                      'baseline_MAE': float(berror.mean()), 'row_coverage': row_coverage,
                      'region_balanced_coverage': region_coverage,
                      'median_halfwidth': float(np.median(width[accepted])),
                      'fraction_radius_at_most_original_width_cap': float(np.mean(width[accepted] <= np.maximum(.2, baseline[accepted])))}
            for key, value in checks.items():
                assert abs(value - expected[key]) < 1e-6, (source, key, value, expected[key])
        audits[source] = {'passed': True, 'actual_saved_decision_rows': len(ids), 'accepted_rows': int(accepted.sum()),
                          'independent_metrics': checks, 'source_result_sha256': sha(result_path),
                          'decision_sha256': sha(decision_path)}
    assert audits
    receipt = {'passed_for_completed_sources': True, 'audited_sources': audits,
        'all_planned_sources_complete_and_audited': len(audits) == len(plan['outer_studies']),
        'scope': 'Actual saved identity, answer/missingness masks and measured-error/coverage arithmetic only. Does not independently reconstruct fitted error critics or calibration factors.',
        'new_point_accuracy_or_useful_confidence_established': False, 'new_study_outcomes_read': False,
        'code_sha256': sha(Path(__file__)), 'plan_sha256': sha(plan_path)}
    (OUT / 'selected_confidence_independent_metric_audit.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt), flush=True)


if __name__ == '__main__':
    main()

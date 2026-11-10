"""Read-only decision-bottleneck audit of the archived confidence policy.

Use only existing development calibration identities and saved predictions.
No thresholds change, model fitting, new outcomes or new validation claims.
"""
import collections
import hashlib
import json
from pathlib import Path
import time
import joblib
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'outputs/sleep_frailty_ml_next'
OLD = DATA / 'validation'
OUT = ROOT / 'outputs/sleep_frailty_ml_confidence_revision'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def read(path):
    return json.loads(path.read_text())


def main():
    start = time.monotonic()
    OUT.mkdir(exist_ok=True)
    settings = read(OLD / 'final_settings_lock.json')
    assessment = read(OLD / 'sequence_confidence_assessment.json')
    audit = read(OLD / 'sequence_confidence_independent_audit.json')
    package_path = OLD / 'sequence_confidence.joblib'
    expected = assessment['confidence_package_sha256']
    assert audit['passed'] and audit['confidence_package_sha256'] == sha(package_path) == expected
    # This is our own already-hash-verified historical sklearn artifact.
    package = joblib.load(package_path)
    policy_path = OLD / 'confidence_policy_lock.json'
    assert package['policy_lock_sha256'] == sha(policy_path)
    policy = read(policy_path)
    assert sha(DATA / 'coverage_training_arrays.npz') == settings['source_data_sha256']
    with np.load(DATA / 'coverage_training_arrays.npz', allow_pickle=False) as z:
        d = {k: z[k] for k in ['split', 'bank_index', 'cell', 'line', 'study', 'group', 'study_names']}
    # Memory-map original arrays; build raw inputs only for calibration rows.
    bank0 = np.load(DATA / 'cell_matched_native_development_features.npy', mmap_mode='r')
    bank1 = np.load(DATA / 'coverage_appended_reader_features.npy', mmap_mode='r')
    phys0 = np.load(DATA / 'native_development_physical_metadata.npy', mmap_mode='r')
    phys1 = np.load(DATA / 'coverage_appended_physical_metadata.npy', mmap_mode='r')
    feature_receipts = [read(DATA / 'cell_matched_feature_receipt.json'), read(DATA / 'coverage_feature_receipt.json')]
    for name, expected_hash in [
        ('cell_matched_native_development_features.npy', feature_receipts[0]['features_sha256']),
        ('native_development_physical_metadata.npy', feature_receipts[0]['physical_metadata_sha256']),
        ('coverage_appended_reader_features.npy', feature_receipts[1]['features_sha256']),
        ('coverage_appended_physical_metadata.npy', feature_receipts[1]['physical_metadata_sha256'])]:
        assert sha(DATA / name) == expected_hash
    all_features, all_base, all_cells, all_cs, all_ls, all_counts, all_ids = [], [], [], [], [], [], []
    source_hashes = {}
    for source in settings['held_studies'][1:]:
        meta, predictions = None, []
        for seed in settings['seeds']:
            path = OLD / 'fits' / f'fixed_{source}_sequence_{seed}_record.json'
            record = read(path)
            assert record['settings_lock_sha256'] == sha(OLD / 'final_settings_lock.json')
            source_hashes[str(path.relative_to(ROOT))] = sha(path)
            p = Path(record['predictions'])
            assert sha(p) == record['predictions_sha256']
            source_hashes[str(p.relative_to(ROOT))] = sha(p)
            with np.load(p, allow_pickle=False) as z:
                mask = d['split'][z['observation_index']] == 2
                item = {k: z[k][mask] for k in z.files}
            if meta is None:
                meta = record
            else:
                for key in ['cell_medians', 'line_medians', 'supported_cells', 'supported_lines', 'training_cell_source_counts', 'training_line_source_counts']:
                    assert record[key] == meta[key]
                assert np.array_equal(item['observation_index'], predictions[0]['observation_index'])
            predictions.append(item)
        ix = predictions[0]['observation_index']
        sid = d['study_names'].tolist().index(source)
        assert np.all(d['study'][ix] == sid)
        means = {k: np.mean([p[k] for p in predictions], axis=0) for k in ['magnitude', 'signed_effect', 'activity_reference', 'activity_alternate']}
        ensemble_sd = np.std([p['magnitude'] for p in predictions], axis=0)
        banks, cell, line = d['bank_index'][ix], d['cell'][ix], d['line'][ix]
        lower = banks < len(bank0)
        bank = np.empty((len(ix), 3, 16), dtype='float32')
        physical = np.empty((len(ix), 4), dtype='float32')
        bank[lower], bank[~lower] = bank0[banks[lower]], bank1[banks[~lower] - len(bank0)]
        physical[lower], physical[~lower] = phys0[banks[lower]], phys1[banks[~lower] - len(phys0)]
        x = np.c_[bank.reshape(len(ix), 48), physical[:, 0] / 230.,
                  physical[:, 1] / np.maximum(physical[:, 0] - 1, 1),
                  np.eye(16, dtype='float32')[(physical[:, 2] * 4 + physical[:, 3]).astype('int64')]].astype('float32')
        cm, lm = np.asarray(meta['cell_medians'])[cell], np.asarray(meta['line_medians'])[line]
        cs, ls = np.asarray(meta['supported_cells'])[cell], np.asarray(meta['supported_lines'])[line]
        counts = np.asarray(meta['training_cell_source_counts'])[cell]
        lcounts = np.asarray(meta['training_line_source_counts'])[line]
        baseline = np.where(ls, lm, cm)
        features = np.column_stack([means['magnitude'], ensemble_sd, abs(means['signed_effect']),
             (means['activity_reference'] + means['activity_alternate']) * .5, cm, baseline,
             counts, lcounts, cs, ls, x, np.eye(len(meta['cell_medians']), dtype='float32')[cell]]).astype('float32')
        assert np.isfinite(features).all()
        all_features.append(features[cs]); all_base.append(baseline[cs]); all_cells.append(cell[cs])
        all_cs.append(cs[cs]); all_ls.append(ls[cs]); all_counts.append(counts[cs]); all_ids.append(ix[cs])
    features, baseline, cell, cs, ls, counts, ids = [np.concatenate(v) for v in [all_features, all_base, all_cells, all_cs, all_ls, all_counts, all_ids]]
    assert len(np.unique(ids)) == len(ids) == assessment['final_scale_rows']
    assert hashlib.sha256(ids.tobytes()).hexdigest() == package['scale_observation_index_sha256']
    critics = package['critics']
    raw = np.maximum(critics['quantile_error'].predict(features), 1e-4)
    estimated = np.maximum(critics['mean_error'].predict(features), 0)
    base = np.maximum(critics['baseline_error'].predict(features), 1e-4)
    scales = [package['by_cell'].get(int(c), package['global_scale']) for c in cell]
    floor = np.array([v['absolute_error_floor'] for v in scales])
    scaled = raw * np.array([v['factor'] for v in scales])
    width = np.maximum(scaled, floor)
    cap = np.maximum(policy['maximum_halfwidth_absolute_floor'], baseline)
    z = np.max(abs((features - package['feature_center']) / package['feature_sd']), axis=1)
    conditions = {'known_cell': cs, 'known_exact_line': ls,
                  'minimum_training_study_count': counts >= policy['minimum_training_studies_same_cell'],
                  'inside_input_support_bound': z <= policy['maximum_standardized_input'],
                  'estimated_relative_gain': (base - estimated) / base >= policy['minimum_estimated_relative_gain'],
                  'estimated_absolute_gain': base - estimated >= policy['minimum_estimated_absolute_gain'],
                  'maximum_interval_halfwidth': width <= cap}
    accepted = np.logical_and.reduce(list(conditions.values()))
    other = np.logical_and.reduce([v for k, v in conditions.items() if k != 'maximum_interval_halfwidth'])
    def summarize(mask):
        return {'rows': int(mask.sum()), 'pretest_policy_accepted': int((mask & accepted).sum()),
            'failures_by_condition_overlap': {k: int((mask & ~v).sum()) for k, v in conditions.items()},
            'passing_all_conditions_except_width': int((mask & other).sum()),
            'unconditional_floor_alone_exceeds_width_cap': int((mask & (floor > cap)).sum()),
            'passing_other_conditions_but_floor_makes_acceptance_impossible': int((mask & other & (floor > cap)).sum()),
            'median_raw_conditional_radius': float(np.median(raw[mask])),
            'median_scaled_conditional_radius': float(np.median(scaled[mask])),
            'median_unconditional_error_floor': float(np.median(floor[mask])),
            'median_permitted_radius': float(np.median(cap[mask]))}
    result = {'passed_identity_and_policy_reconstruction': True,
        'scope': 'Read-only decision-bottleneck diagnosis on archived development calibration rows. No threshold ablation is promoted and no outcome accuracy is calculated.',
        'overall': summarize(np.ones(len(ids), dtype=bool)),
        'by_source': {str(source): summarize(d['study'][ids] == sid) for sid, source in enumerate(d['study_names']) if np.any(d['study'][ids] == sid)},
        'all_calibration_rows_have_exact_original_identity_hash': True,
        'policy_sha256': sha(policy_path), 'confidence_package_sha256': expected,
        'archived_prediction_and_metadata_sha256': source_hashes,
        'raw_effect_labels_or_reserved_outcomes_loaded': False,
        'new_models_fitted': False, 'old_thresholds_or_predictions_changed': False,
        'new_useful_confidence_established': False, 'scientific_goal_complete': False,
        'code_sha256': sha(Path(__file__)), 'elapsed_seconds': time.monotonic() - start}
    (OUT / 'archived_confidence_bottleneck_diagnostic.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['overall', 'by_source', 'elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()

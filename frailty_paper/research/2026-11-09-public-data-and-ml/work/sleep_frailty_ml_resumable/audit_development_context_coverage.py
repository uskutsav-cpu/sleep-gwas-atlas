"""Describe actual source/cell support and target units without changing fits.

This diagnostic uses only the established development corpus, excluding original
calibration rows from target summaries. It does not read the reserved new study,
fit an oracle/test-derived baseline, select a model or alter any current target.
"""
import datetime
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/sleep_frailty_ml_resumable'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()


def main():
    destination = OUT / 'development_context_coverage_audit.json'
    assert not destination.exists()
    path = ROOT / 'outputs/sleep_frailty_ml_next/coverage_training_arrays.npz'
    plan = json.loads((OUT / 'finite_stable_native_plan.json').read_text())
    assert sha(path) == plan['data_sha256']
    keys = ['study', 'cell', 'line', 'context', 'group', 'split', 'effect',
            'activity_reference', 'activity_alternate', 'measurement_standard_error',
            'study_names', 'cell_names', 'cell_line_names', 'context_names']
    with np.load(path, allow_pickle=False) as archive:
        d = {key: archive[key] for key in keys}
    assert np.isfinite(d['effect']).all() and np.all(d['effect'] >= 0)
    active = np.unique(d['context'])
    contexts = d['context_names'][active].tolist()
    assert all(context.endswith('|absolute_log2_effect') for context in contexts)
    assert not any('saturation_coefficient' in context for context in contexts)
    assert len(d['effect']) == 654612
    supported_roles = plan['phase_one_roles'][1:] + plan['phase_two_roles_if_gate_passes']
    roles = {}
    for source in supported_roles:
        sid = d['study_names'].tolist().index(source)
        source_groups = np.unique(d['group'][d['study'] == sid])
        train = (d['split'] == 0) & (d['study'] != sid) & ~np.isin(d['group'], source_groups)
        observed = (d['study'] == sid) & (d['split'] != 2)
        assert not np.intersect1d(d['group'][train], d['group'][observed]).size
        assert not np.any(train & (d['study'] == sid))
        tr, te = np.flatnonzero(train), np.flatnonzero(observed)
        training_hash = hashlib.sha256(tr.tobytes()).hexdigest()
        prediction_hash = hashlib.sha256(te.tobytes()).hexdigest()
        matched_started = []
        for recipe in plan['recipes']:
            metadata_path = OUT / 'development' / f'{recipe}_{source}_{plan["seed"]}_started.json'
            if metadata_path.exists():
                metadata = json.loads(metadata_path.read_text())
                assert metadata['training_index_sha256'] == training_hash
                assert metadata['prediction_index_sha256'] == prediction_hash
                matched_started.append(metadata_path.name)
        cell_support = np.isin(d['cell'][te], np.unique(d['cell'][tr]))
        line_support = np.isin(d['line'][te], np.unique(d['line'][tr]))
        details = []
        for line in np.unique(d['line'][te]):
            selected = observed & (d['line'] == line)
            parent = np.unique(d['cell'][selected])
            assert len(parent) == 1
            cell = int(parent[0])
            same_line = train & (d['line'] == line)
            same_cell = train & (d['cell'] == cell)
            activities = selected & np.isfinite(d['activity_reference']) & np.isfinite(d['activity_alternate'])
            se = selected & np.isfinite(d['measurement_standard_error']) & (d['measurement_standard_error'] > 0)
            details.append({
                'exact_line_name': str(d['cell_line_names'][line]), 'broad_cell_group': str(d['cell_names'][cell]),
                'noncalibration_source_rows': int(selected.sum()),
                'other_source_training_rows_same_exact_line': int(same_line.sum()),
                'other_source_training_rows_same_broad_cell': int(same_cell.sum()),
                'other_training_sources_same_exact_line': d['study_names'][np.unique(d['study'][same_line])].tolist(),
                'other_training_sources_same_broad_cell': d['study_names'][np.unique(d['study'][same_cell])].tolist(),
                'paired_measured_activity_rows': int(activities.sum()),
                'known_positive_measurement_SE_rows': int(se.sum()),
                'source_noncalibration_magnitude_median_for_description_only': float(np.median(d['effect'][selected])),
                'source_noncalibration_magnitude_q90_for_description_only': float(np.quantile(d['effect'][selected], .9)),
                'train_only_line_or_parent_magnitude_median': float(np.median(d['effect'][same_line if same_line.any() else same_cell])) if same_cell.any() else None})
        roles[source] = {
            'training_rows': len(tr), 'noncalibration_source_rows': len(te),
            'exact_line_supported_fraction': float(line_support.mean()),
            'broad_cell_supported_fraction': float(cell_support.mean()),
            'source_absent_and_all_its_region_groups_excluded_from_training': True,
            'training_index_sha256': training_hash, 'prediction_index_sha256': prediction_hash,
            'current_fit_metadata_identities_matched': matched_started, 'cell_line_details': details}
    result = {
        'passed': True, 'time_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'actual_active_observations': len(d['effect']), 'actual_active_contexts': contexts,
        'all_active_primary_targets_have_explicit_log2_effect_unit': True,
        'old_saturation_context_names_retained_but_no_saturation_observations_active': True,
        'same_nominal_unit_does_not_prove_identical_assay_scale_or_measurement_process': True,
        'roles': roles,
        'current_model_cell_context': 'Supported broad-cell and exact-line categorical flags, plus three reader activity outputs in the activity/direction readout. Magnitude also uses the existing 236 DNA/physical features, including projected channels. There is no measured transcriptome representation of a new cell state in this model.',
        'interpretation': 'This describes actual coverage and possible experiment/cell-state transfer limitations. It does not show that adding expression features will help, that source shifts caused the poor performance, or that nominal log2 units make experiments interchangeable.',
        'source_label_summaries_used_as_inference_inputs': False,
        'target_source_median_fitted_as_test_prediction': False,
        'calibration_outcomes_used_for_target_summaries': False,
        'reserved_new_study_outcomes_read': False, 'models_fitted_or_selected': False,
        'current_training_targets_weights_partitions_or_code_changed': False,
        'data_sha256': sha(path), 'code_sha256': sha(Path(__file__))}
    assert sha(path) == plan['data_sha256']
    destination.write_text(json.dumps(result, indent=2) + '\n')
    lines = [
        'Development cell/experiment coverage: checked without changing training', '',
        'All 654,612 active target measurements are explicitly classified as log2 effect magnitudes.',
        'The old saturation-assay context names remain in metadata, but their incompatible-unit rows are absent.',
        'Sharing a unit does not make different experiments or cell maturation states identical.', '',
        'Fraction of omitted-study observations with the same named line in other-source training:']
    for source, row in roles.items():
        lines.append(f'{source}: exact line {row["exact_line_supported_fraction"]:.1%}; broad cell group {row["broad_cell_supported_fraction"]:.1%}')
    lines += ['', 'These are descriptions of the established development data, not new-study validation.',
              'Source outcome summaries are not inputs or fitted test baselines. Calibration outcomes were excluded.',
              'The current fit is unchanged. A missing biological cell-state description is a hypothesis to investigate,',
              'not a demonstrated cause of failure or proof that another model will succeed.']
    (OUT / 'DEVELOPMENT_CONTEXT_COVERAGE.txt').write_text('\n'.join(lines) + '\n')
    print(json.dumps({'passed': True, 'all_active_primary_targets_explicit_log2': True,
                      'support': {name: {key: value[key] for key in ['exact_line_supported_fraction', 'broad_cell_supported_fraction']}
                                  for name, value in roles.items()}, 'models_fitted_or_selected': False}), flush=True)


if __name__ == '__main__':
    main()

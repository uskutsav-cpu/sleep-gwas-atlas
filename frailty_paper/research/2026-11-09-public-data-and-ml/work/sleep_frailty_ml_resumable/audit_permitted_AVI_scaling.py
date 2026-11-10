"""Actual-fold input audit, including held-out perturbation isolation."""
import hashlib
import json
from pathlib import Path
import numpy as np
from permitted_AVI_features import fit_scaling, transform

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/sleep_frailty_ml_resumable/AVI_preparation'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    audit_path = OUT / 'corrected_permitted_AVI_independent_audit.json'
    audit = json.loads(audit_path.read_text())
    assert audit['passed']
    array = OUT / 'corrected_permitted_AVI_physical_bank.npy'
    assert sha(array) == audit['array_sha256']
    bank = np.load(array, allow_pickle=False)
    data = ROOT / 'outputs/sleep_frailty_ml_next/coverage_training_arrays.npz'
    plan = json.loads((ROOT / 'outputs/sleep_frailty_ml_resumable/finite_stable_native_plan.json').read_text())
    assert sha(data) == plan['data_sha256']
    with np.load(data, allow_pickle=False) as z:
        d = {key: z[key] for key in ['bank_index', 'split', 'study', 'study_names', 'group']}
    roles = plan['phase_one_roles'] + plan['phase_two_roles_if_gate_passes']
    assert len(roles) == len(set(roles)) == 8
    records = {}
    for role in roles:
        train = d['split'] == 0
        if role != 'none':
            sid = d['study_names'].tolist().index(role)
            groups = np.unique(d['group'][d['study'] == sid])
            train &= (d['study'] != sid) & ~np.isin(d['group'], groups)
        ids = np.unique(d['bank_index'][train])
        settings = fit_scaling(bank, ids)
        reference = bank[ids, :2].astype('float64')
        assert np.isfinite(reference).all()
        assert np.allclose(settings['median'], np.quantile(reference, .5, axis=0), atol=0, rtol=0)
        assert np.allclose(settings['mean'], reference.sum(axis=0) / len(reference), atol=1e-12, rtol=0)
        variance = ((reference - reference.mean(axis=0)) ** 2).sum(axis=0) / len(reference)
        assert np.allclose(settings['sd'], np.sqrt(variance), atol=1e-12, rtol=0)
        # Counterfactual engineering probe, never written as data or used to
        # train a model: change ONLY nontraining inputs, then refit scaling.
        copy = bank.copy()
        nontraining = np.setdiff1d(np.unique(d['bank_index']), ids)
        assert len(nontraining)
        copy[nontraining, :2] = [1e6, -1e6]
        assert fit_scaling(copy, ids) == settings
        del copy
        probe = bank[np.r_[ids[:3], 279177]].copy()
        with_scores = transform(probe, settings, True)
        without_scores = transform(probe, settings, False)
        assert np.array_equal(with_scores[:, 2], without_scores[:, 2])
        assert np.all(without_scores[:, :2] == 0)
        assert probe[-1, 0] == 0 and probe[-1, 2] == 1 and np.isfinite(with_scores[-1]).all()
        absent = transform(np.array([[np.nan, np.nan, 0]], dtype='float32'), settings, True)
        assert np.isfinite(absent).all() and np.all(absent[0, :2] == 0)
        assert absent[0, 2] == -settings['training_availability_mean']
        original = ROOT / 'outputs/sleep_frailty_ml_resumable/development' / f'frozen_reader_{role}_42_started.json'
        index_sha = hashlib.sha256(np.flatnonzero(train).tobytes()).hexdigest()
        if original.exists():
            assert json.loads(original.read_text())['training_index_sha256'] == index_sha
        records[role] = {'scaling': settings, 'training_observation_index_sha256': index_sha,
                         'held_out_input_perturbations_changed_scaling': False,
                         'both_arms_have_identical_availability_features': True,
                         'valid_zero_is_finite_and_known': True, 'missing_inputs_keep_an_explicit_flag': True}
    result = {'passed': True, 'roles': records,
        'comparison_design': 'Same 3-column head shape in both arms. Only normalized raw/PHRED score values differ; availability is shared. Initialize all three extra input weight columns to zero and preserve common parameters, sampling, updates and randomness.',
        'scoreless_control_description': 'DNA predictor with the same availability metadata; no composite score values.',
        'unknown_input_policy': 'Raw array retains missing values. The model sees neutral standardized score features and an explicit availability difference; no zero biological-effect label is imputed.',
        'new_models_fitted_or_promoted': False, 'reserved_outcomes_or_effect_labels_read': False,
        'array_sha256': sha(array), 'independent_input_audit_sha256': sha(audit_path),
        'code_sha256': {p.name: sha(p) for p in [Path(__file__), Path(__file__).with_name('permitted_AVI_features.py')]},
        'remaining': 'Implement and audit the matched head/training experiment once its sequence-model recipe and development gates are fixed. Freeze both arms and final-study inputs/predictions before the new outcome comparison.'}
    (OUT / 'permitted_AVI_train_only_scaling_audit.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'passed': True, 'roles': len(records), 'new_models_fitted': False, 'held_out_inputs_affect_training_scaling': False}), flush=True)


if __name__ == '__main__':
    main()

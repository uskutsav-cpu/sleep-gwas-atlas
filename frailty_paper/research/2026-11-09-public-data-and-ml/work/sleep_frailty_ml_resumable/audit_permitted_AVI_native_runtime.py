"""Replay complete historical DNA constructs through the actual new model.

This is a discarded numerical audit, not a candidate fit or new-study test.
"""
import copy
import json
from pathlib import Path
import resource
import shutil

from resumable_common import ROOT, assert_versions_stable
from native_data import PackedDNA
from bounded_network import NativeAdaptationModel, torch, np
from permitted_AVI_network import PermittedAVIModel
from permitted_AVI_head import expand_head_for_scores
from safetensors.torch import load_file
from readonly_native_inputs import sha
from audit_permitted_AVI_head import numpy_forward, check_native_forward_derivation

OUT = ROOT / 'outputs/sleep_frailty_ml_resumable'
AVI = OUT / 'AVI_preparation'


def main():
    destination = AVI / 'permitted_AVI_native_runtime_audit.json'
    assert not destination.exists()
    assert shutil.disk_usage(ROOT).free > 128 * 1024**2
    assert_versions_stable()
    check_native_forward_derivation()
    torch.set_num_threads(1)
    case_receipt = json.loads((AVI / 'native_runtime_diagnostic_cases_receipt.json').read_text())
    path = AVI / 'native_runtime_diagnostic_cases.npz'
    assert sha(path) == case_receipt['diagnostic_data_sha256']
    with np.load(path, allow_pickle=False) as archive:
        cases = {key: archive[key] for key in archive.files}
    n = len(cases['length'])
    packed = PackedDNA.__new__(PackedDNA)
    packed.length, packed.position, packed.alt = cases['length'], cases['position'], cases['alt']
    packed.offset = np.r_[0, np.cumsum(packed.length.astype('int64'))]
    packed.tokens = cases['reference_tokens']
    assert len(packed.tokens) == packed.offset[-1]
    assert sorted(set(map(int, packed.length))) == [150, 199, 200, 270]
    d = dict(bank_index=np.arange(n, dtype='int64'), cell=cases['cell'], line=cases['line'])
    physical = np.zeros((n, 66), dtype='float32')
    physical[:, 48:50] = cases['physical']
    name = 'frozen_reader_none_42'
    directory = OUT / 'development'
    meta = json.loads((directory / (name + '_started.json')).read_text())
    assert sha(directory / (name + '_started.json')) == case_receipt['original_model_metadata_sha256']
    original_receipt = json.loads((directory / (name + '_epoch3_receipt.json')).read_text())
    checkpoint = directory / (name + '_epoch3.safetensors')
    assert sha(checkpoint) == original_receipt['checkpoint_sha256']
    head_state = {key.removeprefix('head.'): value for key, value in load_file(str(checkpoint)).items()}
    args = (np.arange(n), d, packed, cases['bank'], physical, cases['channels'],
            np.asarray(meta['supported_cells']), np.asarray(meta['supported_lines']))
    original = NativeAdaptationModel(meta['cell_medians'], meta['line_medians'], meta['line_parent'],
                                     meta['normalization'], 'frozen_reader').double().eval()
    original.head.load_state_dict(head_state, strict=True)
    with torch.no_grad():
        expected = {key: value.clone() for key, value in original(*args).items()}
    fitted_head = copy.deepcopy(original.head)
    del original
    with np.load(directory / (name + '_epoch3_predictions.npz'), allow_pickle=False) as saved:
        ids = saved['observation_index']
        assert sha(directory / (name + '_epoch3_predictions.npz')) == original_receipt['prediction_sha256']
        positions = np.searchsorted(ids, cases['observation_ids'][:8])
        assert np.array_equal(ids[positions], cases['observation_ids'][:8])
        saved_errors = {key: float(abs(expected[key].numpy()[:8] - saved[key][positions]).max()) for key in expected}
    assert max(saved_errors.values()) < 1e-4
    scaling = json.loads((AVI / 'permitted_AVI_train_only_scaling_audit.json').read_text())['roles']['none']['scaling']
    equivalence, numpy_errors, gradients = {}, {}, []
    for use_scores in [False, True]:
        model = PermittedAVIModel(meta['cell_medians'], meta['line_medians'], meta['line_parent'],
                                  meta['normalization'], 'frozen_reader', scaling, use_scores).double().eval()
        model.head = expand_head_for_scores(fitted_head)
        captured = []
        hook = model.head.register_forward_pre_hook(lambda module, values: captured.append(values))
        actual = model(*args, cases['avi'])
        hook.remove()
        error = max(float((actual[key] - expected[key]).abs().max().detach()) for key in actual)
        assert error < 1e-12
        independent = numpy_forward(model.head, captured[0])
        ne = {key: float(abs(actual[key].detach().numpy() - independent[key]).max()) for key in actual}
        assert max(ne.values()) < 1e-12
        arm = 'scores' if use_scores else 'scoreless'
        equivalence[arm], numpy_errors[arm] = error, ne
        if use_scores:
            target = torch.from_numpy(cases['effect'].astype('float64'))
            def objective():
                return (model(*args, cases['avi'])['magnitude'] - target).abs().mean()
            objective().backward()
            parameter = model.head.hidden.weight
            values = parameter.grad[:, 236:238]
            flat = int(values.abs().argmax())
            row, column = flat // 2, 236 + flat % 2
            analytic = float(parameter.grad[row, column])
            assert abs(analytic) > 1e-10
            epsilon = 1e-6
            with torch.no_grad():
                before = float(parameter[row, column])
                parameter[row, column] = before + epsilon
                plus = float(objective())
                parameter[row, column] = before - epsilon
                minus = float(objective())
                parameter[row, column] = before
            finite = (plus - minus) / (2 * epsilon)
            relative = abs(finite - analytic) / max(abs(finite), abs(analytic), 1e-12)
            assert relative < 1e-5
            gradients.append({'parameter': 'head.hidden.weight', 'row': row, 'column': column,
                              'relative_gradient_error': relative})
            # A fixed discarded perturbation verifies that the new connection
            # can alter a prediction without breaking allele-swap symmetry.
            with torch.no_grad():
                parameter[:, 236:238] = .002
                forward = model(*args, cases['avi'])
                reverse = model(*args, cases['avi'], swap=True)
                assert torch.allclose(forward['magnitude'], reverse['magnitude'], atol=1e-12, rtol=0)
                assert torch.allclose(forward['signed_effect'], -reverse['signed_effect'], atol=1e-12, rtol=0)
                response = float((forward['magnitude'] - expected['magnitude']).abs().max())
                assert response > 1e-10
                parameter[:, 236:238].zero_()
        del model, captured, actual
    assert_versions_stable()
    result = {'passed': True, 'full_native_forward_runtime_replay_complete': True,
              'actual_historical_observations_checked': n, 'actual_DNA_lengths': sorted(set(map(int, packed.length))),
              'full_actual_REF_and_ALT_sequences_used': True, 'invented_or_padded_bases': False,
              'original_fitted_head_checkpoint_sha256': original_receipt['checkpoint_sha256'],
              'expanded_models_vs_original_max_errors': equivalence,
              'independent_NumPy_head_errors_on_actual_DNA_inputs': numpy_errors,
              'CPU_float64_vs_saved_MPS_float32_max_errors': saved_errors, 'saved_prediction_absolute_tolerance': 1e-4,
              'new_input_finite_gradient_checks_on_real_DNA_and_historical_effects': gradients,
              'fixed_discarded_score_weight_perturbation_changes_prediction': response,
              'nonzero_score_connection_keeps_allele_swap_size_invariance_and_signed_antisymmetry': True,
              'actual_zero_score_and_reversed_reference_role_cases_included': True,
              'no_optimizer_steps_or_scientific_model_selection': True,
              'original_saved_weights_and_scientific_results_unchanged': True,
              'reserved_new_study_outcomes_read': False, 'predictive_benefit_established': False,
              'peak_RSS_bytes_on_macOS': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              'case_receipt_sha256': sha(AVI / 'native_runtime_diagnostic_cases_receipt.json'),
              'code_sha256': {p.name: sha(p) for p in [Path(__file__), Path(__file__).with_name('permitted_AVI_network.py'), Path(__file__).with_name('permitted_AVI_head.py')]},
              'remaining': 'Fix the candidate and the matched scientific training/selection protocol; fit both arms, then freeze final-study predictions before outcome evaluation. Numerical checks do not prove a prediction improvement.'}
    destination.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in ['passed', 'actual_historical_observations_checked', 'actual_DNA_lengths', 'expanded_models_vs_original_max_errors', 'peak_RSS_bytes_on_macOS', 'predictive_benefit_established']}), flush=True)


if __name__ == '__main__':
    main()

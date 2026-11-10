"""Actual fitted-head/gradient checks, without loading a second DNA reader.

The feature vectors in this engineering audit are fixed diagnostic inputs, not
experimental training examples. No scientific model is fitted or selected.
"""
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil

import numpy as np
import torch
from safetensors.torch import load_file

from permitted_AVI_features import transform
from permitted_AVI_head import expand_head_for_scores

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/sleep_frailty_ml_resumable'
WORK = Path(__file__).parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def forward_tree(path, classname):
    tree = ast.parse(path.read_text())
    cls = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == classname)
    return next(node for node in cls.body if isinstance(node, ast.FunctionDef) and node.name == 'forward')


def check_native_forward_derivation():
    old = forward_tree(ROOT / 'work/sleep_frailty_ml_bounded/bounded_network.py', 'NativeAdaptationModel')
    new = forward_tree(WORK / 'permitted_AVI_network.py', 'PermittedAVIModel')
    assert [arg.arg for arg in new.args.args].count('avi_bank') == 1
    new.args.args = [arg for arg in new.args.args if arg.arg != 'avi_bank']
    added = [ast.parse(line).body[0] for line in [
        "score_features=torch.from_numpy(transform(avi_bank[b],self.score_scaling,self.use_scores)).to(device=device,dtype=dtype)",
        "features=torch.cat([features,score_features],1)"]]
    wanted = [ast.dump(node) for node in added]
    removed = [node for node in new.body if ast.dump(node) in wanted]
    assert len(removed) == 2
    new.body = [node for node in new.body if ast.dump(node) not in wanted]
    assert ast.dump(old) == ast.dump(new), 'An unrelated native forward computation changed.'


def numpy_forward(head, args):
    state = {key: value.detach().numpy() for key, value in head.state_dict().items()}
    x, rr, aa, cell, cf, lf = [v.detach().numpy() for v in args]
    weight = state['cell_activity_weight'][cell] + lf @ state['line_activity_weight']
    bias = state['cell_activity_bias'][cell] + lf @ state['line_activity_bias']
    reference = (rr * weight).sum(1) + bias
    alternate = (aa * weight).sum(1) + bias
    signed = alternate - reference
    inputs = np.column_stack([x, abs(signed) * 10., (reference + alternate) * .5, cf, lf])
    def silu(values):
        return values / (1 + np.exp(np.clip(-values, -600, 600)))
    hidden = silu(inputs @ state['hidden.weight'].T + state['hidden.bias'])
    hidden = silu(hidden @ state['hidden2.weight'].T + state['hidden2.bias'])
    magnitude = np.logaddexp(0., state['cell_bias'][cell] + lf @ state['line_offset'] +
                             hidden @ state['output.weight'][0] + state['output.bias'][0])
    return dict(magnitude=magnitude, signed_effect=signed,
                activity_reference=reference, activity_alternate=alternate)


def main():
    assert shutil.disk_usage(ROOT).free > 256 * 1024**2
    torch.set_num_threads(1)
    check_native_forward_derivation()
    name = 'frozen_reader_none_42'
    directory = OUT / 'development'
    meta = json.loads((directory / (name + '_started.json')).read_text())
    receipt = json.loads((directory / (name + '_epoch3_receipt.json')).read_text())
    checkpoint = directory / (name + '_epoch3.safetensors')
    assert sha(checkpoint) == receipt['checkpoint_sha256']
    head_source = ROOT / 'work/sleep_frailty_ml_next/cell_matched_effect_network.py'
    spec = importlib.util.spec_from_file_location('audited_original_cell_head', head_source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    base = module.CellMatchedEffectHead(meta['cell_medians'], meta['line_medians'], meta['line_parent'], 236)
    state = load_file(str(checkpoint))
    assert all(key.startswith('head.') for key in state)
    base.load_state_dict({key.removeprefix('head.'): value for key, value in state.items()}, strict=True)
    base.double().eval()
    original = {key: value.clone() for key, value in base.state_dict().items()}
    rng = torch.get_rng_state().clone()
    control = expand_head_for_scores(base).eval()
    scored = expand_head_for_scores(base).eval()
    assert torch.equal(rng, torch.get_rng_state())
    for key in control.state_dict():
        assert torch.equal(control.state_dict()[key], scored.state_dict()[key])
        if key != 'hidden.weight':
            assert torch.equal(control.state_dict()[key], original[key])
    for head in [control, scored]:
        assert torch.equal(head.hidden.weight[:, :236], base.hidden.weight[:, :236])
        assert torch.count_nonzero(head.hidden.weight[:, 236:239]) == 0
        assert torch.equal(head.hidden.weight[:, 239:], base.hidden.weight[:, 236:])
    preparation = OUT / 'AVI_preparation'
    bank = np.load(preparation / 'corrected_permitted_AVI_physical_bank.npy', mmap_mode='r')
    scaling_receipt = json.loads((preparation / 'permitted_AVI_train_only_scaling_audit.json').read_text())
    assert scaling_receipt['passed'] and sha(preparation / 'corrected_permitted_AVI_physical_bank.npy') == scaling_receipt['array_sha256']
    physical = np.r_[np.flatnonzero(bank[:, 2] == 1)[:16],
                     np.flatnonzero((bank[:, 2] == 1) & (bank[:, 0] == 0))]
    values = np.array(bank[physical], copy=True)
    values = np.vstack([values, [np.nan, np.nan, 0]])
    scaling = scaling_receipt['roles']['none']['scaling']
    cfeatures = transform(values, scaling, False)
    sfeatures = transform(values, scaling, True)
    assert np.array_equal(cfeatures[:, 2], sfeatures[:, 2])
    assert not cfeatures[:, :2].any() and np.isfinite(sfeatures).all()
    n = len(values)
    features = np.sin(np.arange(n * 236).reshape(n, 236) / 117.) * .3
    rr = np.sin(np.arange(n * 3).reshape(n, 3) / 9.)
    aa = rr + np.cos(np.arange(n * 3).reshape(n, 3) / 11.) * .02
    parent = np.asarray(meta['line_parent'])
    line = np.arange(n) % len(parent)
    cell = parent[line]
    cf = np.eye(len(meta['cell_medians']))[cell]
    lf = np.eye(len(parent))[line]
    tensor = lambda value: torch.from_numpy(np.asarray(value))
    suffix = [tensor(rr), tensor(aa), tensor(cell), tensor(cf), tensor(lf)]
    original_args = [tensor(features), *suffix]
    cargs = [tensor(np.column_stack([features, cfeatures])), *suffix]
    sargs = [tensor(np.column_stack([features, sfeatures])), *suffix]
    expected = base(*original_args)
    assert float(expected['magnitude'].detach().std()) > 1e-8
    initial_error = 0.
    for head, args in [(control, cargs), (scored, sargs)]:
        prediction = head(*args)
        for key in expected:
            error = float((prediction[key] - expected[key]).abs().max().detach())
            initial_error = max(initial_error, error)
            assert error < 1e-12
    # The same dropout mask must also give the same initial predictions.
    control.train(); scored.train()
    state_rng = torch.get_rng_state().clone()
    cp = control(*cargs)
    torch.set_rng_state(state_rng)
    sp = scored(*sargs)
    assert max(float((cp[key] - sp[key]).abs().max().detach()) for key in cp) < 1e-12
    control.eval(); scored.eval()
    # Diagnostic targets test the real backward implementation, not usefulness.
    target = tensor(.12 + .03 * np.sin(np.arange(n)))
    def loss(head, args):
        return (head(*args)['magnitude'] - target).abs().mean()
    closs, sloss = loss(control, cargs), loss(scored, sargs)
    closs.backward(); sloss.backward()
    common_gradient_error = 0.
    for key, param in control.named_parameters():
        other = dict(scored.named_parameters())[key]
        cg, sg = param.grad, other.grad
        assert (cg is None) == (sg is None)
        if cg is None:
            continue
        assert torch.isfinite(cg).all() and torch.isfinite(sg).all()
        if key == 'hidden.weight':
            cg = torch.cat([cg[:, :236], cg[:, 239:]], 1)
            sg = torch.cat([sg[:, :236], sg[:, 239:]], 1)
        error = float((cg - sg).abs().max())
        common_gradient_error = max(common_gradient_error, error)
        assert error < 1e-12
    assert not control.hidden.weight.grad[:, 236:238].any()
    region = scored.hidden.weight.grad[:, 236:238]
    flat = int(region.abs().argmax())
    row, column = flat // 2, 236 + flat % 2
    analytic = float(scored.hidden.weight.grad[row, column])
    assert abs(analytic) > 1e-10
    epsilon = 1e-6
    with torch.no_grad():
        original_value = scored.hidden.weight[row, column].item()
        scored.hidden.weight[row, column] = original_value + epsilon
        plus = float(loss(scored, sargs))
        scored.hidden.weight[row, column] = original_value - epsilon
        minus = float(loss(scored, sargs))
        scored.hidden.weight[row, column] = original_value
    finite = (plus - minus) / (2 * epsilon)
    relative = abs(finite - analytic) / max(abs(finite), abs(analytic), 1e-12)
    assert relative < 1e-5
    for head in [control, scored]:
        torch.optim.AdamW(head.parameters(), lr=3e-4, weight_decay=.001).step()
    assert torch.count_nonzero(scored.hidden.weight[:, 236:238]) > 0
    assert torch.count_nonzero(control.hidden.weight[:, 236:238]) == 0
    numpy_errors = {}
    for name, head, args in [('scoreless', control, cargs), ('scored', scored, sargs)]:
        actual, independent = head(*args), numpy_forward(head, args)
        errors = {key: float(np.max(abs(actual[key].detach().numpy() - independent[key]))) for key in actual}
        assert max(errors.values()) < 1e-12
        numpy_errors[name] = errors
        swapped = head(args[0], args[2], args[1], *args[3:])
        assert torch.allclose(actual['magnitude'], swapped['magnitude'], atol=1e-12, rtol=0)
        assert torch.allclose(actual['signed_effect'], -swapped['signed_effect'], atol=1e-12, rtol=0)
    scored_response = float((scored(*sargs)['magnitude'] - scored(*cargs)['magnitude']).abs().max().detach())
    assert scored_response > 1e-10
    assert all(torch.equal(original[key], value) for key, value in base.state_dict().items())
    assert sha(checkpoint) == receipt['checkpoint_sha256']
    result = {
        'passed': True, 'native_forward_unchanged_except_two_score_concatenation_statements': True,
        'original_head_checkpoint_sha256': receipt['checkpoint_sha256'],
        'identical_between_arm_initial_states_and_CPU_RNG_preserved': True,
        'all_original_parameters_preserved': True, 'same_dropout_mask_preserves_initial_predictions': True,
        'initial_prediction_max_error': initial_error, 'common_initial_gradient_max_error': common_gradient_error,
        'new_score_weight_finite_gradient_relative_error': relative,
        'score_weights_receive_gradients_and_change_on_diagnostic_step': True,
        'scoreless_raw_and_PHRED_weight_columns_remain_zero': True,
        'score_input_prediction_response_after_diagnostic_step': scored_response,
        'independent_NumPy_forward_max_errors': numpy_errors,
        'allele_swap_magnitude_invariance_and_signed_antisymmetry': True,
        'head_parameters_per_arm': sum(p.numel() for p in scored.parameters()),
        'added_parameter_scalars_per_arm': 3 * scored.hidden.out_features,
        'known_zero_score_examples_checked': int(np.sum((values[:, 2] == 1) & (values[:, 0] == 0))),
        'availability_inputs_identical_between_arms_including_unknown_case': True,
        'diagnostic_examples_not_actual_experimental_training_examples': True,
        'scientific_models_trained_or_selected': False, 'reserved_outcomes_read': False,
        'full_native_forward_runtime_replay_complete': False,
        'scope': 'Actual expanded fitted-head forward/backward and one discarded diagnostic optimizer step; AST proof that the native forward differs only by the stated score-feature concatenation. No predictive-benefit evidence.',
        'remaining': 'Full native runtime replay, fixed candidate/training/selection protocol, matched scientific fits, and frozen untouched-study benefit evaluation.',
        'code_sha256': {p.name: sha(p) for p in [Path(__file__), WORK / 'permitted_AVI_head.py', WORK / 'permitted_AVI_network.py', WORK / 'permitted_AVI_features.py', head_source]}}
    destination = preparation / 'permitted_AVI_expanded_head_audit.json'
    assert not destination.exists(), 'Preserve prior audit evidence rather than overwriting it.'
    destination.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in ['passed', 'initial_prediction_max_error', 'common_initial_gradient_max_error', 'new_score_weight_finite_gradient_relative_error', 'head_parameters_per_arm', 'scientific_models_trained_or_selected']}), flush=True)


if __name__ == '__main__':
    main()

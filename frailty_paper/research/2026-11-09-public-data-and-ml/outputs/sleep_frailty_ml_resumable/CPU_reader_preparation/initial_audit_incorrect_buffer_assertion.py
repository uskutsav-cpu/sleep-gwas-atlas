"""Small actual-model numerical/device/recovery check; no scientific fitting."""
import hashlib
import json
from pathlib import Path
import resource
import subprocess

from resumable_common import ROOT
from bounded_network import NativeAdaptationModel, torch, np
from native_data import PackedDNA
from cell_matched_effect_network import CellMatchedEffectHead
from safetensors.torch import load_file
from cpu_reader_mps_head import CPUReaderMPSHead
from resumption import load_recovery
from readonly_native_inputs import sha

OUT = ROOT / 'outputs/sleep_frailty_ml_resumable'


def parameter_hash(model):
    digest = hashlib.sha256()
    for name, value in sorted(model.state_dict().items()):
        digest.update(name.encode()); digest.update(value.detach().cpu().numpy().tobytes())
    return digest.hexdigest()


def main():
    destination = OUT / 'CPU_reader_preparation/cpu_reader_numerical_and_restore_audit.json'
    assert not destination.exists()
    running = subprocess.run(['pgrep', '-fl', 'train_resumable_development|run_resumable_comparison'],
                             capture_output=True, text=True)
    assert running.returncode == 1, 'Do not overlap this diagnostic with the old scientific worker.'
    assert torch.backends.mps.is_available()
    assert CPUReaderMPSHead.forward is NativeAdaptationModel.forward
    assert CPUReaderMPSHead.train is NativeAdaptationModel.train
    cases_path = OUT / 'AVI_preparation/native_runtime_diagnostic_cases.npz'
    case_receipt = json.loads(cases_path.with_name('native_runtime_diagnostic_cases_receipt.json').read_text())
    assert sha(cases_path) == case_receipt['diagnostic_data_sha256']
    with np.load(cases_path, allow_pickle=False) as archive:
        cases = {key: archive[key] for key in archive.files}
    n = len(cases['length'])
    packed = PackedDNA.__new__(PackedDNA)
    packed.length, packed.position, packed.alt = cases['length'], cases['position'], cases['alt']
    packed.offset = np.r_[0, np.cumsum(packed.length.astype('int64'))]
    packed.tokens = cases['reference_tokens']
    physical = np.zeros((n, 66), dtype='float32'); physical[:, 48:50] = cases['physical']
    meta = json.loads((OUT / 'development/adapt_last_stage_none_42_started.json').read_text())
    model_args = [meta['cell_medians'], meta['line_medians'], meta['line_parent'], meta['normalization'], 'adapt_last_stage']
    args = [np.arange(n), {'bank_index': np.arange(n), 'cell': cases['cell'], 'line': cases['line']},
            packed, cases['bank'], physical, cases['channels'], np.asarray(meta['supported_cells']), np.asarray(meta['supported_lines'])]
    reference = NativeAdaptationModel(*model_args).to('cpu').eval()
    candidate = CPUReaderMPSHead(*model_args).to('mps').eval()
    assert parameter_hash(reference) == parameter_hash(candidate) == meta['initial_parameter_sha256']
    checkpoint = OUT / 'development/adapt_last_stage_none_42_epoch1.safetensors'
    checkpoint_receipt = json.loads(checkpoint.with_name(checkpoint.stem + '_receipt.json').read_text())
    assert sha(checkpoint) == checkpoint_receipt['checkpoint_sha256']
    state = load_file(str(checkpoint))
    for model in [reference, candidate]:
        missing, unexpected = model.load_state_dict(state, strict=False)
        assert not unexpected
        assert all(not key.startswith('head.') and not key.startswith('reader.main.3.') for key in missing)
    assert parameter_hash(reference) == parameter_hash(candidate)
    actual, expected = candidate(*args), reference(*args)
    errors = {key: float((actual[key].detach().cpu() - expected[key].detach().cpu()).abs().max()) for key in actual}
    assert max(errors.values()) < 1e-4, errors
    target = torch.from_numpy(cases['effect'])
    loss = (actual['magnitude'] - target.to('mps')).abs().mean() + candidate.anchor_loss()
    cpu_loss = (expected['magnitude'] - target).abs().mean() + reference.anchor_loss()
    assert abs(float(loss.detach().cpu()) - float(cpu_loss.detach())) < 1e-5
    loss.backward(); cpu_loss.backward()
    gradients = []
    for name, param in candidate.named_parameters():
        other = dict(reference.named_parameters())[name]
        assert (param.grad is None) == (other.grad is None), name
        if param.grad is None:
            continue
        got, wanted = param.grad.detach().cpu().double(), other.grad.detach().cpu().double()
        absolute = float((got - wanted).abs().max())
        relative = float(torch.linalg.vector_norm(got - wanted) / torch.linalg.vector_norm(wanted).clamp_min(1e-8))
        assert absolute < 2e-5 and (relative < .005 or float(wanted.abs().max()) < 1e-6), (name, absolute, relative)
        gradients.append({'parameter': name, 'maximum_absolute_error': absolute, 'relative_L2_error': relative})
    candidate.train()
    standalone = CellMatchedEffectHead(meta['cell_medians'], meta['line_medians'], meta['line_parent'], 236).to('mps').train()
    standalone.load_state_dict(candidate.head.state_dict(), strict=True)
    captured = []
    hook = candidate.head.register_forward_pre_hook(lambda module, values: captured.append(tuple(value.detach() for value in values)))
    torch.mps.synchronize(); before = torch.mps.get_rng_state().clone()
    with torch.no_grad():
        first = candidate(*args)
        after = torch.mps.get_rng_state().clone()
        torch.mps.set_rng_state(before)
        second = standalone(*captured[-1])
        after_second = torch.mps.get_rng_state().clone()
    hook.remove()
    assert torch.equal(after, after_second)
    dropout_error = max(float((first[key] - second[key]).abs().max().cpu()) for key in first)
    assert dropout_error < 1e-6
    # Restore the real interrupted frozen-head state without changing its file
    # or taking a training step. The old helper verifies every parameter, Adam
    # value, and CPU/MPS RNG state after loading into the mixed-device model.
    current = OUT / 'development/frozen_reader_Kosicki2025_42_recovery.pt'
    interruption = json.loads((OUT / 'verified_storage_interruption_Kosicki2025_20261109.json').read_text())
    assert sha(current) == interruption['archived_recovery_sha256']
    saved = json.loads((current.with_name('frozen_reader_Kosicki2025_42_started.json')).read_text())
    frozen = CPUReaderMPSHead(saved['cell_medians'], saved['line_medians'], saved['line_parent'], saved['normalization'], 'frozen_reader').to('mps')
    assert parameter_hash(frozen) == saved['initial_parameter_sha256']
    optimizer = torch.optim.AdamW([{'params': list(frozen.head.parameters()), 'lr': .0003}], weight_decay=.001)
    keys = ['recipe', 'held_study', 'seed', 'epochs', 'training_index_sha256', 'prediction_index_sha256',
            'initial_parameter_sha256', 'plan_sha256', 'code_sha256']
    identity = {key: saved[key] for key in keys}
    identity['recovery_helper_sha256'] = sha(Path(__file__).with_name('resumption.py'))
    restored = load_recovery(current, frozen, optimizer, 'mps', identity)
    assert restored['restore_verified_exactly'] and restored['updates'] == 2600
    assert sha(current) == interruption['archived_recovery_sha256']
    result = {'passed': True, 'original_forward_and_train_methods_inherited_unchanged': True,
              'initial_and_fitted_parameter_state_bytes_identical': True,
              'actual_DNA_examples': n, 'actual_DNA_lengths': np.unique(cases['length']).tolist(),
              'reader_and_feature_devices': sorted({str(p.device) for p in candidate.reader.parameters()}),
              'head_parameter_devices': sorted({str(p.device) for p in candidate.head.parameters()}),
              'CPU_reference_vs_CPU_reader_MPS_head_prediction_errors': errors,
              'full_primary_plus_anchor_gradient_checks': gradients,
              'head_dropout_random_state_consumption_preserved_exactly': True,
              'standalone_head_vs_routed_head_dropout_prediction_max_error': dropout_error,
              'actual_old_checkpoint_parameters_Adam_CPU_and_MPS_rng_restored_exactly': True,
              'old_checkpoint_update': restored['updates'], 'old_checkpoint_file_changed': False,
              'training_updates_or_scientific_fits_performed': False,
              'full_batch_training_resource_or_continuation_preflight_complete': False,
              'macOS_peak_RSS_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              'MPS_allocated_bytes_after_small_audit': torch.mps.current_allocated_memory(),
              'MPS_driver_allocated_bytes_after_small_audit': torch.mps.driver_allocated_memory(),
              'scope': 'Small actual-model arithmetic, gradients, dropout-stream and existing-state restoration only. Rounding across devices is tolerated at stated bounds, not claimed bitwise identical. A pinned execution adapter and full-batch sustained preflight are still required before scientific continuation.',
              'reserved_new_study_outcomes_read': False, 'model_promoted': False,
              'code_sha256': {path.name: sha(path) for path in [Path(__file__), Path(__file__).with_name('cpu_reader_mps_head.py')]}}
    destination.parent.mkdir(exist_ok=True)
    destination.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in ['passed', 'CPU_reference_vs_CPU_reader_MPS_head_prediction_errors',
                     'head_dropout_random_state_consumption_preserved_exactly',
                     'actual_old_checkpoint_parameters_Adam_CPU_and_MPS_rng_restored_exactly',
                     'full_batch_training_resource_or_continuation_preflight_complete', 'macOS_peak_RSS_bytes']}), flush=True)


if __name__ == '__main__':
    main()

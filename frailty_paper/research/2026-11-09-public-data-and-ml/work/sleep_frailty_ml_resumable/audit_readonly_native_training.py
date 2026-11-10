"""Compare actual four-update states with the original MPS reference runs."""
import argparse
import datetime
import json
from pathlib import Path

from audit_recovery_state import compare, float32_replay, torch
from resumable_common import ROOT
from readonly_native_inputs import sha

OUT = ROOT / 'outputs/sleep_frailty_ml_resumable'
RUNS = OUT / 'recovery_audit_runs'
SUFFIX = 'development/adapt_last_stage_none_42_recovery.pt'


def compare_states(reference_path, other_path):
    a = torch.load(reference_path, map_location='cpu', weights_only=True)
    b = torch.load(other_path, map_location='cpu', weights_only=True)
    differences = []
    parameters = float32_replay(a['parameters'], b['parameters'], 'parameters', differences)
    optimizer = float32_replay(a['optimizer'], b['optimizer'], 'optimizer', differences)
    for key in ['CPU_rng', 'MPS_rng']:
        compare(a[key], b[key], key)
    pa, pb = json.loads(a['progress_json']), json.loads(b['progress_json'])
    for key in ['updates', 'epoch', 'phase', 'schedule_sha256', 'NumPy_rng_json', 'curve']:
        compare(pa[key], pb[key], key)
    for key in ['objectives', 'anchors']:
        float32_replay(torch.tensor(pa[key], dtype=torch.float32),
                       torch.tensor(pb[key], dtype=torch.float32), key, differences)
    assert pa['updates'] == pb['updates'] == 4
    return {'reference_sha256': sha(reference_path), 'compared_sha256': sha(other_path),
            'parameter_scalars': parameters, 'optimizer_scalars': optimizer,
            'CPU_and_MPS_RNG_exact': True, 'NumPy_and_batch_schedule_exact': True,
            'numeric_bound': 'Existing eight-float32-ULP four-update bound; initial restore is checked exactly.',
            'differences': differences}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--require-continuation', action='store_true')
    args = parser.parse_args()
    direct = RUNS / 'readonly_native_adapt_direct4'
    source = json.loads((direct / 'effective_runtime_composition.json').read_text())
    for name, digest in source['code_sha256'].items():
        path = (ROOT / 'work/sleep_frailty_ml_bounded' / name if name == 'bounded_network.py'
                else Path(__file__).with_name(name))
        assert sha(path) == digest, name
    original = RUNS / 'reference_four_steps' / SUFFIX
    results = {'readonly_vs_original': compare_states(original, direct / SUFFIX),
               'original_repeat_background_variation': compare_states(original, RUNS / 'repeat_reference_four_steps' / SUFFIX)}
    keys = ['training_index_sha256', 'prediction_index_sha256', 'initial_parameter_sha256',
            'normalization', 'context_transformations', 'new_original_Siraj_signed_training_rows',
            'new_original_neuronal_signed_training_rows']
    original_meta = json.loads(original.with_name('adapt_last_stage_none_42_started.json').read_text())
    direct_meta = json.loads((direct / SUFFIX).with_name('adapt_last_stage_none_42_started.json').read_text())
    for key in keys:
        assert original_meta[key] == direct_meta[key], key
    continuation = RUNS / 'readonly_native_adapt_continued4'
    if args.require_continuation:
        results['resumed_vs_direct'] = compare_states(direct / SUFFIX, continuation / SUFFIX)
        history_path = continuation / 'development/adapt_last_stage_none_42_resume_history.json'
        history = json.loads(history_path.read_text())
        assert history[-1]['updates'] == 2 and history[-1]['model_optimizer_and_CPU_MPS_random_state_restored_bit_for_bit']
        results['exact_initial_restoration'] = {'passed': True, 'history_sha256': sha(history_path)}
    result = {'passed': True, 'time_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'scope': 'Four actual full-batch MPS updates with original model and all losses; not sustained-resource proof or biological validation.',
              'continuation_checked': args.require_continuation, 'matched_metadata_fields': keys,
              'comparisons': results, 'scientific_weights_changed': False,
              'code_sha256': {p.name: sha(p) for p in [Path(__file__), Path(__file__).with_name('audit_recovery_state.py')]},
              'runtime_composition_sha256': sha(direct / 'effective_runtime_composition.json')}
    destination = OUT / ('readonly_native_recovery_audit.json' if args.require_continuation else 'readonly_native_direct_four_update_audit.json')
    destination.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'passed': True, 'continuation_checked': args.require_continuation,
                      'report': str(destination), 'max_observed_ULP': max((row['maximum_float32_ULP'] for r in results.values() for row in r.get('differences', [])), default=0)}), flush=True)


if __name__ == '__main__':
    main()

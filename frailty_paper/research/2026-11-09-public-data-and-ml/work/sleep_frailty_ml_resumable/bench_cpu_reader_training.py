"""Bounded engineering runs of the unchanged trainer with explicit device routing.

No scientific continuation is authorized by this benchmark. Its original-source
metadata is supplemented by a mandatory runtime-composition receipt. All runs
live in the old trainer's separate recovery_audit_runs directory and must stop
at an explicit positive update limit. Their learned weights are discarded.
"""
import argparse
import datetime
import importlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys

os.environ.setdefault('PYTORCH_MPS_HIGH_WATERMARK_RATIO', '0.4')
os.environ.setdefault('PYTORCH_MPS_LOW_WATERMARK_RATIO', '0.3')

from resumable_common import ROOT
from cpu_reader_mps_head import CPUReaderMPSHead, torch
from readonly_native_inputs import load_readonly_inputs, sha

OUT = ROOT / 'outputs/sleep_frailty_ml_resumable'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--recipe', choices=['frozen_reader', 'adapt_last_stage'], required=True)
    parser.add_argument('--audit-tag', required=True)
    parser.add_argument('--stop-after-updates', type=int, required=True)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    assert args.audit_tag.startswith('cpu_reader_') and args.audit_tag.replace('_', '').isalnum()
    assert 0 < args.stop_after_updates <= 1000
    p = subprocess.run(['pgrep', '-fl', 'train_resumable_development|run_resumable_comparison'],
                       capture_output=True, text=True)
    assert p.returncode == 1, 'Do not overlap this resource diagnostic with a scientific native worker.'
    proof_path = OUT / 'CPU_reader_preparation/cpu_reader_numerical_and_restore_audit.json'
    proof = json.loads(proof_path.read_text())
    assert proof['passed']
    for name, digest in proof['code_sha256'].items():
        assert sha(Path(__file__).with_name(name)) == digest, name
    original = importlib.import_module('train_resumable_development')
    original_path = Path(original.__file__)
    engineering = OUT / 'recovery_audit_runs' / args.audit_tag
    engineering.mkdir(exist_ok=True, parents=True)
    composition = {
        'time_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'engineering_only': True,
        'PID': os.getpid(), 'recipe': args.recipe, 'held_study': 'none',
        'runtime_class': 'CPUReaderMPSHead', 'reader_and_features_device': 'cpu', 'head_and_data_loss_device': 'mps',
        'data_loader': 'Audited read-only batch indexing; original values and training normalization unchanged.',
        'original_trainer_metadata_is_not_the_full_runtime_source_inventory': True,
        'learned_audit_weights_authorized_for_scientific_initialization_or_selection': False,
        'original_checkpoint_or_scientific_fit_changed': False,
        'code_sha256': {path.name: sha(path) for path in [Path(__file__), original_path,
                            Path(__file__).with_name('cpu_reader_mps_head.py'),
                            Path(__file__).with_name('readonly_native_inputs.py'),
                            Path(__file__).with_name('resumption.py')]},
        'numerical_proof_sha256': sha(proof_path)}
    receipt = engineering / 'effective_runtime_composition.json'
    if args.resume:
        old = json.loads(receipt.read_text())
        assert old['code_sha256'] == composition['code_sha256']
        assert old['numerical_proof_sha256'] == composition['numerical_proof_sha256']
    else:
        assert not receipt.exists()
        receipt.write_text(json.dumps(composition, indent=2) + '\n')
    original.NativeAdaptationModel = CPUReaderMPSHead
    original.load_data = lambda: load_readonly_inputs()[:5]
    real_save = original.save_recovery
    profile_path = engineering / ('resource_profile_resumed.json' if args.resume else 'resource_profile.json')
    assert not profile_path.exists()
    profile = []
    def save_with_resource_profile(path, model, optimizer, device, identity, progress):
        real_save(path, model, optimizer, device, identity, progress)
        profile.append({'updates': progress['updates'], 'phase': progress['phase'],
                        'maximum_RSS_bytes_on_macOS': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                        'MPS_current_allocated_bytes': torch.mps.current_allocated_memory(),
                        'MPS_driver_allocated_bytes': torch.mps.driver_allocated_memory(),
                        'actual_reader_parameter_devices': sorted({str(p.device) for p in model.reader.parameters()}),
                        'actual_head_parameter_devices': sorted({str(p.device) for p in model.head.parameters()})})
        profile_path.write_text(json.dumps(profile, indent=2) + '\n')
    original.save_recovery = save_with_resource_profile
    sys.argv = [str(original_path), '--recipe', args.recipe, '--held-study', 'none',
                '--audit-tag', args.audit_tag, '--audit-stop-after-updates', str(args.stop_after_updates)] + (['--resume'] if args.resume else [])
    original.main()
    name = f'{args.recipe}_none_42'
    stopped = json.loads((engineering / 'development' / (name + '_audit_stop.json')).read_text())
    assert stopped['updates'] == args.stop_after_updates and not stopped['used_to_initialize_scientific_fit']
    result = {'bounded_engineering_segment_completed': True, 'actual_updates': stopped['updates'],
              'MPS_head_CPU_reader_full_batch_training_executed': True,
              'scientific_model_or_goal_success_established': False,
              'effective_runtime_composition_sha256': sha(receipt), 'resource_profile_sha256': sha(profile_path)}
    result_path = engineering / ('resumed_segment_result.json' if args.resume else 'initial_segment_result.json')
    result_path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()

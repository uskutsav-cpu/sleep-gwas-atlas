"""Bounded original-MPS training audit with only the input allocation replaced.

The model, loss, device, batches and optimizer are the original implementation.
Audit weights stay in recovery_audit_runs and cannot initialize scientific fits.
The original trainer's source inventory is explicitly supplemented here.
"""
import argparse
import datetime
import importlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

os.environ.setdefault('PYTORCH_MPS_HIGH_WATERMARK_RATIO', '0.4')
os.environ.setdefault('PYTORCH_MPS_LOW_WATERMARK_RATIO', '0.3')

from resumable_common import ROOT
from readonly_native_inputs import load_readonly_inputs, sha

OUT = ROOT / 'outputs/sleep_frailty_ml_resumable'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--recipe', choices=['frozen_reader', 'adapt_last_stage'], required=True)
    parser.add_argument('--audit-tag', required=True)
    parser.add_argument('--stop-after-updates', type=int, required=True)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    assert args.audit_tag.startswith('readonly_native_') and args.audit_tag.replace('_', '').isalnum()
    assert 0 < args.stop_after_updates <= 1000
    processes = subprocess.run(['ps', '-axo', 'pid,command'], check=True, capture_output=True, text=True).stdout
    for line in processes.splitlines()[1:]:
        pid, command = line.strip().split(None, 1)
        if int(pid) != os.getpid() and 'python' in command:
            assert not any(str(Path(__file__).with_name(name)) in command for name in
                           ['train_resumable_development.py', 'run_resumable_comparison.py', Path(__file__).name]), line
    assert shutil.disk_usage(ROOT).free >= 1024**3, 'Need at least 1 GiB before a bounded resource audit.'
    proof_path = OUT / 'AVI_preparation/readonly_native_input_equivalence_audit.json'
    proof = json.loads(proof_path.read_text())
    assert proof['passed'] and proof['full_training_mean_and_SD_bitwise_equal_to_original_fit']
    for name, digest in proof['code_sha256'].items():
        assert sha(Path(__file__).with_name(name)) == digest, name
    original = importlib.import_module('train_resumable_development')
    original_path = Path(original.__file__)
    torch = original.torch
    engineering = OUT / 'recovery_audit_runs' / args.audit_tag
    engineering.mkdir(exist_ok=True, parents=True)
    composition = {
        'time_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'engineering_only': True,
        'PID': os.getpid(), 'recipe': args.recipe, 'held_study': 'none',
        'runtime_class': original.NativeAdaptationModel.__name__, 'reader_and_head_device': 'mps',
        'data_loader': 'Read-only batch indexing; original physical and reader values unchanged.',
        'original_trainer_metadata_is_not_the_full_runtime_source_inventory': True,
        'learned_audit_weights_authorized_for_scientific_initialization_or_selection': False,
        'original_checkpoint_or_scientific_fit_changed': False,
        'code_sha256': {path.name: sha(path) for path in [Path(__file__), original_path,
                            Path(__file__).with_name('readonly_native_inputs.py'),
                            Path(__file__).with_name('resumption.py'),
                            Path(inspect_module_file(original.NativeAdaptationModel))]},
        'input_equivalence_proof_sha256': sha(proof_path)}
    receipt = engineering / 'effective_runtime_composition.json'
    if args.resume:
        old = json.loads(receipt.read_text())
        for key in ['code_sha256', 'input_equivalence_proof_sha256']:
            assert old[key] == composition[key], key
    else:
        assert not receipt.exists()
        receipt.write_text(json.dumps(composition, indent=2) + '\n')
    original.load_data = lambda: load_readonly_inputs()[:5]
    real_save = original.save_recovery
    profile_path = engineering / ('resource_profile_resumed.json' if args.resume else 'resource_profile.json')
    assert not profile_path.exists()
    profile = []
    started = time.monotonic()

    def save_with_resource_profile(path, model, optimizer, device, identity, progress):
        real_save(path, model, optimizer, device, identity, progress)
        profile.append({'updates': progress['updates'], 'phase': progress['phase'],
                        'segment_elapsed_seconds': time.monotonic() - started,
                        'maximum_RSS_bytes_on_macOS': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                        'free_system_disk_bytes': shutil.disk_usage(ROOT).free,
                        'MPS_current_allocated_bytes': torch.mps.current_allocated_memory(),
                        'MPS_driver_allocated_bytes': torch.mps.driver_allocated_memory(),
                        'actual_reader_parameter_devices': sorted({str(p.device) for p in model.reader.parameters()}),
                        'actual_head_parameter_devices': sorted({str(p.device) for p in model.head.parameters()})})
        profile_path.write_text(json.dumps(profile, indent=2) + '\n')
        assert profile[-1]['free_system_disk_bytes'] > 300*1024**2, 'Resource audit stops before disk reserve is exhausted.'

    original.save_recovery = save_with_resource_profile
    sys.argv = [str(original_path), '--recipe', args.recipe, '--held-study', 'none',
                '--audit-tag', args.audit_tag, '--audit-stop-after-updates', str(args.stop_after_updates)] + (['--resume'] if args.resume else [])
    original.main()
    name = f'{args.recipe}_none_42'
    stopped = json.loads((engineering / 'development' / (name + '_audit_stop.json')).read_text())
    assert stopped['updates'] == args.stop_after_updates and not stopped['used_to_initialize_scientific_fit']
    result = {'bounded_engineering_segment_completed': True, 'actual_updates': stopped['updates'],
              'original_MPS_model_and_full_batch_training_executed': True,
              'only_data_allocation_changed': True,
              'scientific_model_or_goal_success_established': False,
              'segment_elapsed_seconds': time.monotonic() - started,
              'effective_runtime_composition_sha256': sha(receipt), 'resource_profile_sha256': sha(profile_path)}
    result_path = engineering / ('resumed_segment_result.json' if args.resume else 'initial_segment_result.json')
    result_path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


def inspect_module_file(cls):
    return sys.modules[cls.__module__].__file__


if __name__ == '__main__':
    main()

"""Admit only a completed actual 1000-update original-MPS allocation audit."""
import datetime
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'outputs/sleep_frailty_ml_resumable'
RUN = OUT/'recovery_audit_runs/readonly_native_adapt_resource1000'


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024**2),b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    composition_path = RUN/'effective_runtime_composition.json'
    composition = json.loads(composition_path.read_text())
    actual = subprocess.run(['ps','-p',str(composition['PID']),'-o','command='],capture_output=True,text=True)
    assert 'bench_readonly_native_training.py' not in actual.stdout, 'Wait for the existing resource process to exit.'
    result_path = RUN/'initial_segment_result.json'
    result = json.loads(result_path.read_text())
    assert result['bounded_engineering_segment_completed'] and result['actual_updates'] == 1000
    assert result['original_MPS_model_and_full_batch_training_executed'] and result['only_data_allocation_changed']
    assert result['effective_runtime_composition_sha256'] == sha(composition_path)
    profile_path = RUN/'resource_profile.json'
    assert result['resource_profile_sha256'] == sha(profile_path)
    profile = json.loads(profile_path.read_text())
    assert [r['updates'] for r in profile] == list(range(0,1001,100))
    assert all(r['actual_reader_parameter_devices'] == ['mps:0'] and r['actual_head_parameter_devices'] == ['mps:0'] for r in profile)
    minimum_free = min(r['free_system_disk_bytes'] for r in profile)
    peak_rss = max(r['maximum_RSS_bytes_on_macOS'] for r in profile)
    late_growth = profile[-1]['maximum_RSS_bytes_on_macOS'] - profile[5]['maximum_RSS_bytes_on_macOS']
    assert minimum_free > 300*1024**2, 'Reserve breached at a recorded checkpoint.'
    assert profile[-1]['free_system_disk_bytes'] >= 1024**3, 'Not enough headroom for a controlled continuation.'
    assert peak_rss < 3*1024**3 and late_growth < 256*1024**2
    for name,digest in composition['code_sha256'].items():
        base = ROOT/'work/sleep_frailty_ml_bounded' if name == 'bounded_network.py' else Path(__file__).parent
        assert sha(base/name) == digest, name
    four_path = OUT/'readonly_native_recovery_audit.json'
    four = json.loads(four_path.read_text())
    assert four['passed'] and four['continuation_checked']
    input_path = OUT/'AVI_preparation/readonly_native_input_equivalence_audit.json'
    assert composition['input_equivalence_proof_sha256'] == sha(input_path)
    resource_result = {'passed':True,'time_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                       'actual_training_updates':1000,'actual_resource_process_terminal':True,
                       'original_model_and_device':True,'maximum_RSS_bytes_on_macOS':peak_rss,
                       'minimum_free_disk_at_recorded_checkpoints':minimum_free,
                       'final_free_disk_bytes':profile[-1]['free_system_disk_bytes'],
                       'late_500_update_peak_RSS_growth_bytes':late_growth,
                       'largest_recorded_MPS_driver_allocation_bytes':max(r['MPS_driver_allocated_bytes'] for r in profile),
                       'GPU_and_free_disk_samples_are_checkpoint_observations_not_continuous_peak_measurements':True,
                       'scope':'1000 actual full-batch adaptation updates. Original full-validation preflight remains inherited; this is not a new full-validation run or evidence of predictive improvement.',
                       'audit_weights_used_to_initialize_scientific_fits':False,
                       'evidence_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [composition_path,result_path,profile_path,four_path,input_path]},
                       'code_sha256':{str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))}}
    path = OUT/'readonly_native_sustained_resource_audit.json'
    assert not path.exists()
    path.write_text(json.dumps(resource_result,indent=2)+'\n')
    print(json.dumps({key:resource_result[key] for key in ['passed','actual_training_updates','maximum_RSS_bytes_on_macOS','minimum_free_disk_at_recorded_checkpoints','late_500_update_peak_RSS_growth_bytes']}))


if __name__ == '__main__':
    main()

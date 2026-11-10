"""Version only resource admission/receipt paths; preserve the actual fitter.

Source preparation is safe while V1 trains. --admit requires every scientific
worker to be terminal, no pending fit, and a completed observed resource trace.
This is not a new neural numerical audit or a scientific performance gate.
"""
import argparse
import ast
import datetime
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).parent
OUT = ROOT/'outputs/sleep_frailty_ml_resumable'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024**2),b''):
            h.update(block)
    return h.hexdigest()


def replacements(source, changes):
    for old,new in changes:
        assert source.count(old) == 1, old
        source = source.replace(old,new,1)
    return source


def prepare():
    original_revision = json.loads((OUT/'readonly_native_execution_revision.json').read_text())
    runtime = WORK/'readonly_native_runtime.py'
    driver = WORK/'readonly/run_resumable_comparison.py'
    for path in [runtime,driver]:
        assert sha(path) == original_revision['execution_code_sha256'][str(path.relative_to(ROOT))]
    runtime_changes = [
        ("path = DEST/'readonly_native_execution_revision.json'", "path = DEST/'readonly_native_execution_revision_768.json'"),
        ("assert shutil.disk_usage(ROOT).free >= 1024**3, 'Need 1 GiB headroom for this execution revision.'",
         "assert shutil.disk_usage(ROOT).free >= 768*1024**2, 'Need 768 MiB headroom for the admitted resource-only execution revision.'"),
        ("DEST/'readonly_runtime_recovery_audit.json'", "DEST/'readonly_runtime_768_admission_audit.json'"),
        ("Path(__file__).parent/'readonly/run_resumable_comparison.py'", "Path(__file__).parent/'readonly_768/run_resumable_comparison.py'"),
        ("DEST/'readonly_native_recovery_audit.json', DEST/'readonly_native_sustained_resource_audit.json']",
         "DEST/'readonly_native_recovery_audit.json', DEST/'readonly_native_sustained_resource_audit.json',\n              DEST/'resource_admission_768_source_derivation.json']"),
        ("receipt_path = directory/(name+'_effective_runtime.json')", "receipt_path = directory/(name+'_effective_runtime_768.json')"),
        ("names = {'train_resumable_development.py', 'readonly_native_runtime.py',",
         "names = {'train_resumable_development.py', 'readonly_native_runtime.py', 'readonly_native_runtime_768.py',")]
    driver_changes = [
        ('from readonly_native_runtime import execution_revision','from readonly_native_runtime_768 import execution_revision'),
        ("script = str(PARENT/'readonly_native_runtime.py')", "script = str(PARENT/'readonly_native_runtime_768.py')"),
        ("DEST/'readonly_runtime_recovery_audit.json'", "DEST/'readonly_runtime_768_admission_audit.json'"),
        ("DEST/'readonly_finite_driver_execution.json'", "DEST/'readonly_768_finite_driver_execution.json'")]
    outputs = [(runtime,WORK/'readonly_native_runtime_768.py',runtime_changes),
               (driver,WORK/'readonly_768/run_resumable_comparison.py',driver_changes)]
    audits = []
    for original,target,changes in outputs:
        old = original.read_text()
        new = replacements(old,changes)
        ast.parse(new,filename=str(target))
        restored = new
        for before,after in reversed(changes):
            assert restored.count(after) == 1
            restored = restored.replace(after,before,1)
        assert restored == old
        target.parent.mkdir(exist_ok=True)
        if target.exists():
            assert target.read_text() == new, target
        else:
            target.write_text(new)
        audits.append({'original':str(original.relative_to(ROOT)),'original_sha256':sha(original),
                       'derived':str(target.relative_to(ROOT)),'derived_sha256':sha(target),
                       'declared_exact_replacements':changes,'entire_original_source_recovered_by_reversal':True})
    result = {'passed':True,'time_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'scope':'Mechanical source derivation for a resource-only admission change. Not a new neural numerical or prediction audit.',
              'old_launch_reserve_bytes':1024**3,'new_launch_reserve_bytes':768*1024**2,
              'committed_checkpoint_live_reserve_bytes_unchanged':300*1024**2,
              'source_derivations':audits,'original_native_trainer_sha256':sha(WORK/'train_resumable_development.py'),
              'model_loss_labels_sampler_optimizer_seed_epochs_or_scientific_gates_changed':False,
              'current_V1_execution_or_completed_files_changed':False,
              'new_scientific_training_admitted':False,'builder_sha256':sha(Path(__file__))}
    path = OUT/'resource_admission_768_source_derivation.json'
    if path.exists():
        prior = json.loads(path.read_text())
        for key in result:
            if key != 'time_UTC':
                assert prior[key] == result[key], key
    else:
        path.write_text(json.dumps(result,indent=2)+'\n')
    return path


def admit(derivation_path):
    processes = subprocess.run(['ps','-axo','pid,command'],check=True,capture_output=True,text=True).stdout
    worker_names = {'readonly_native_runtime.py','readonly_native_runtime_768.py',
                    'train_resumable_development.py','bench_readonly_native_training.py','bench_cpu_reader_training.py'}
    for row in processes.splitlines()[1:]:
        fields = row.strip().split(None,1)
        if len(fields) != 2 or int(fields[0]) == os.getpid():
            continue
        arguments = shlex.split(fields[1])
        if arguments and Path(arguments[0]).name.startswith('python'):
            assert not any(Path(arg).name in worker_names for arg in arguments[1:]), 'Wait for the actual scientific worker to exit.'
    records = list((OUT/'development').glob('*_record.json'))
    assert len(records) >= 6, 'Finish and preserve the current matched Kosicki comparison first.'
    for path in (OUT/'development').glob('*_started.json'):
        assert path.with_name(path.name.replace('_started.json','_record.json')).exists(), 'No partial fit may migrate through this resource-only admission.'
    trace_path = OUT/'development/adapt_last_stage_Kosicki2025_42_readonly_resource_history.json'
    trace = json.loads(trace_path.read_text())
    assert trace[-1]['phase'] == 'epoch_complete'
    minimum = min(row['free_system_disk_bytes'] for row in trace)
    assert minimum > 300*1024**2
    assert max(row['maximum_RSS_bytes_on_macOS'] for row in trace) < 3*1024**3
    assert shutil.disk_usage(ROOT).free >= 768*1024**2
    sys.path.insert(0,str(WORK))
    from readonly_native_runtime_768 import execution_revision
    revision,revision_hash = execution_revision()
    assert not revision['allowed_original_checkpoint_migrations']
    old_numerics = OUT/'readonly_runtime_recovery_audit.json'
    assert json.loads(old_numerics.read_text())['passed']
    proof = {'passed':True,'time_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
             'execution_revision_sha256':revision_hash,'admitted_launch_reserve_bytes':768*1024**2,
             'live_checkpoint_reserve_bytes_unchanged':300*1024**2,
             'exact_original_runtime_recovered_by_declared_source_reversal':True,
             'new_numerical_training_audit_claimed':False,
             'inherited_actual_numerical_audit_sha256':sha(old_numerics),
             'source_derivation_sha256':sha(derivation_path),
             'completed_adapted_fit_trace_sha256':sha(trace_path),
             'observed_minimum_free_disk_in_completed_adapted_fit':minimum,
             'all_prior_completed_fits_preserved':True,'partial_fit_migration_authorized':False,
             'scope':'Only prelaunch disk admission and versioned receipts change. The actual model, data, optimizer and scientific decision code are inherited unchanged. This is not a new performance result.'}
    path = OUT/'readonly_runtime_768_admission_audit.json'
    assert not path.exists()
    path.write_text(json.dumps(proof,indent=2)+'\n')
    print(json.dumps({'resource_only_revision_admitted':True,'launch_reserve_MiB':768,
                      'checkpoint_reserve_MiB':300,'new_model_numerical_audit_claimed':False}),flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--admit',action='store_true')
    args = parser.parse_args()
    derivation = prepare()
    if args.admit:
        admit(derivation)
    else:
        print(json.dumps({'resource_only_sources_prepared':True,'scientific_admission_complete':False,
                          'must_wait_for_current_worker_terminal':True}),flush=True)


if __name__ == '__main__':
    main()

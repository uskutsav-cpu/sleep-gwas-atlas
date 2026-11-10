"""Attested allocation-only continuation of the original scientific trainer.

No architecture, device, labels, loss, sampling or scientific decision changes.
The inherited source inventory is supplemented, never presented as complete.
New recovery states bind the execution revision into their checked identity.
"""
import argparse
import datetime
import importlib
import json
import os
from pathlib import Path
import resource
import shlex
import shutil
import subprocess
import sys

os.environ.setdefault('PYTORCH_MPS_HIGH_WATERMARK_RATIO', '0.4')
os.environ.setdefault('PYTORCH_MPS_LOW_WATERMARK_RATIO', '0.3')
from resumable_common import ROOT, DEST, assert_versions_stable
from readonly_native_inputs import load_readonly_inputs, sha


def execution_revision():
    assert_versions_stable()
    code = [Path(__file__), Path(__file__).with_name('readonly_native_inputs.py'),
            Path(__file__).with_name('train_resumable_development.py'),
            Path(__file__).with_name('resumption.py'), Path(__file__).with_name('resumable_common.py'),
            Path(__file__).with_name('run_resumable_comparison.py'),
            Path(__file__).with_name('process_control.py'),
            Path(__file__).parent/'readonly/run_resumable_comparison.py']
    proofs = [DEST/'AVI_preparation/readonly_native_input_equivalence_audit.json',
              DEST/'readonly_native_recovery_audit.json', DEST/'readonly_native_sustained_resource_audit.json']
    for path in proofs:
        assert json.loads(path.read_text())['passed'], path
    path = DEST/'readonly_native_execution_revision.json'
    code_hashes = {str(p.relative_to(ROOT)): sha(p) for p in code}
    proof_hashes = {str(p.relative_to(ROOT)): sha(p) for p in proofs}
    if path.exists():
        revision = json.loads(path.read_text())
        assert revision['execution_code_sha256'] == code_hashes
        assert revision['admission_proofs_sha256'] == proof_hashes
    else:
        migration = {}
        for started in (DEST/'development').glob('*_started.json'):
            name = started.name.removesuffix('_started.json')
            if not started.with_name(name+'_record.json').exists():
                recovery = started.with_name(name+'_recovery.pt')
                assert recovery.exists(), name
                migration[name] = {'checkpoint_sha256': sha(recovery), 'started_metadata_sha256': sha(started)}
        revision = {'time_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    'original_scientific_plan_sha256': sha(DEST/'finite_stable_native_plan.json'),
                    'execution_code_sha256': code_hashes, 'admission_proofs_sha256': proof_hashes,
                    'allowed_original_checkpoint_migrations': migration,
                    'change': 'Only replace whole-bank feature allocation with numerically audited read-only batch indexing.',
                    'model_device_targets_objectives_schedule_and_gates_changed': False,
                    'engineering_weights_used_for_scientific_initialization': False}
        path.write_text(json.dumps(revision, indent=2)+'\n')
    assert revision['original_scientific_plan_sha256'] == sha(DEST/'finite_stable_native_plan.json')
    return revision, sha(path)


def assert_no_other_worker():
    rows = subprocess.run(['ps','-axo','pid,command'],check=True,capture_output=True,text=True).stdout
    names = {'train_resumable_development.py', 'readonly_native_runtime.py', 'bench_readonly_native_training.py',
             'bench_cpu_reader_training.py', 'train_permitted_AVI_development.py'}
    for row in rows.splitlines()[1:]:
        fields = row.strip().split(None,1)
        if len(fields) != 2 or int(fields[0]) == os.getpid():
            continue
        arguments = shlex.split(fields[1])
        if arguments and Path(arguments[0]).name.startswith('python'):
            assert not any(Path(arg).name in names for arg in arguments[1:]), ('Concurrent training is forbidden',row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--recipe',choices=['frozen_reader','adapt_last_stage'],required=True)
    parser.add_argument('--held-study',required=True)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--audit-tag')
    parser.add_argument('--audit-stop-after-updates',type=int)
    args = parser.parse_args()
    assert_no_other_worker()
    assert shutil.disk_usage(ROOT).free >= 1024**3, 'Need 1 GiB headroom for this execution revision.'
    if args.audit_tag:
        assert args.audit_tag.startswith('readonly_runtime_') and args.audit_stop_after_updates
        assert 0 < args.audit_stop_after_updates <= 4
    else:
        assert args.audit_stop_after_updates is None
    revision, revision_hash = execution_revision()
    if not args.audit_tag:
        proof = json.loads((DEST/'readonly_runtime_recovery_audit.json').read_text())
        assert proof['passed'] and proof['execution_revision_sha256'] == revision_hash
    original = importlib.import_module('train_resumable_development')
    torch = original.torch
    directory = DEST/'recovery_audit_runs'/args.audit_tag/'development' if args.audit_tag else DEST/'development'
    directory.mkdir(parents=True,exist_ok=True)
    name = f'{args.recipe}_{args.held_study}_42'
    receipt_path = directory/(name+'_effective_runtime.json')
    receipt = {'execution_revision_sha256': revision_hash,
               'actual_execution_overrides_sha256': revision['execution_code_sha256'],
               'inherited_original_trainer_code_inventory_is_partial_without_this_addendum': True,
               'original_model_class_and_forward_unchanged': True, 'reader_and_head_device': 'mps',
               'audit_only': bool(args.audit_tag)}
    if receipt_path.exists():
        assert json.loads(receipt_path.read_text()) == receipt
    else:
        receipt_path.write_text(json.dumps(receipt,indent=2)+'\n')
    original.load_data = lambda: load_readonly_inputs()[:5]
    real_save, real_load, real_write = original.save_recovery, original.load_recovery, original.write_json
    resource_path = directory/(name+'_readonly_resource_history.json')
    resources = json.loads(resource_path.read_text()) if resource_path.exists() else []

    def extended(identity):
        return {**identity,'execution_revision_sha256':revision_hash}

    def save(path, model, optimizer, device, identity, progress):
        real_save(path,model,optimizer,device,extended(identity),progress)
        resources.append({'PID':os.getpid(),'time_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                          'updates':progress['updates'],'phase':progress['phase'],
                          'maximum_RSS_bytes_on_macOS':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                          'MPS_driver_allocated_bytes':torch.mps.driver_allocated_memory(),
                          'free_system_disk_bytes':shutil.disk_usage(ROOT).free})
        real_write(resource_path,resources)
        assert resources[-1]['free_system_disk_bytes'] > 300*1024**2, 'Stop with committed recovery state before exhausting disk.'

    def load(path, model, optimizer, device, identity):
        payload = torch.load(path,map_location='cpu',weights_only=True)
        saved_identity = json.loads(payload['identity_json'])
        del payload
        if 'execution_revision_sha256' in saved_identity:
            assert saved_identity == extended(identity)
            return real_load(path,model,optimizer,device,extended(identity))
        assert not args.audit_tag and args.resume, 'Only the pinned interrupted scientific checkpoint may migrate.'
        allowed = revision['allowed_original_checkpoint_migrations'][name]
        assert sha(path) == allowed['checkpoint_sha256']
        assert sha(directory/(name+'_started.json')) == allowed['started_metadata_sha256']
        assert saved_identity == identity
        progress = real_load(path,model,optimizer,device,identity)
        real_write(directory/(name+'_readonly_migration.json'),
                   {'original_checkpoint_sha256':allowed['checkpoint_sha256'],
                    'execution_revision_sha256':revision_hash,'PID':os.getpid(),
                    'restored_updates':progress['updates'],'restored_epoch':progress['epoch'],
                    'model_optimizer_and_CPU_MPS_RNG_exactly_restored':progress['restore_verified_exactly'],
                    'original_batch_schedule_must_be_verified_before_next_update':True})
        return progress

    def write(path, value):
        if Path(path).name in {name+'_started.json',name+'_record.json'}:
            value = {**value,'execution_revision_sha256':revision_hash,
                     'effective_runtime_addendum_sha256':sha(receipt_path)}
        return real_write(path,value)

    original.save_recovery, original.load_recovery, original.write_json = save, load, write
    sys.argv = [str(original.__file__),'--recipe',args.recipe,'--held-study',args.held_study]
    if args.resume:
        sys.argv.append('--resume')
    if args.audit_tag:
        sys.argv += ['--audit-tag',args.audit_tag,'--audit-stop-after-updates',str(args.audit_stop_after_updates)]
    original.main()


if __name__ == '__main__':
    main()

"""Inspect the completed execution-bound recovery state after its worker exits."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'outputs/sleep_frailty_ml_resumable'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda:handle.read(1024**2),b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--recipe',required=True)
    parser.add_argument('--source',required=True)
    args = parser.parse_args()
    processes = subprocess.run(['ps','-axo','command'],check=True,capture_output=True,text=True).stdout
    worker_names = ['readonly_native_runtime.py --recipe ', 'readonly_native_runtime_768.py --recipe ']
    assert not any(any(name in row for name in worker_names) and 'python' in row for row in processes.splitlines()), 'Do not overlap checkpoint inspection and training.'
    import torch
    from safetensors.torch import load_file
    name = f'{args.recipe}_{args.source}_42'
    base = OUT/'development'
    metadata = json.loads((base/(name+'_started.json')).read_text())
    record_path = base/(name+'_record.json')
    record = json.loads(record_path.read_text())
    revision_candidates = [path for path in OUT.glob('readonly_native_execution_revision*.json')
                           if sha(path) == record['execution_revision_sha256']]
    assert revision_candidates, 'No preserved execution revision matches the completed fit.'
    revision_path = sorted(revision_candidates)[0]
    revision_hash = sha(revision_path)
    revision = json.loads(revision_path.read_text())
    for path,digest in revision['execution_code_sha256'].items():
        assert sha(ROOT/path) == digest,path
    recovery = base/(name+'_recovery.pt')
    state = torch.load(recovery,map_location='cpu',weights_only=True)
    identity = json.loads(state['identity_json'])
    assert identity['execution_revision_sha256'] == record['execution_revision_sha256'] == revision_hash
    for key,value in identity.items():
        if key == 'execution_revision_sha256':
            continue
        if key == 'recovery_helper_sha256':
            assert sha(Path(__file__).with_name('resumption.py')) == value
        else:
            assert metadata[key] == value,key
    progress = json.loads(state['progress_json'])
    assert progress['phase'] == 'epoch_complete'
    assert progress['epoch'] == record['completed_epochs'] == len(progress['curve']) == 3
    assert progress['updates'] == record['actual_backprop_updates']
    assert progress['schedule_sha256'] == record['actual_batch_schedule_sha256']
    final_receipt_path = base/(name+'_epoch3_receipt.json')
    final_receipt = json.loads(final_receipt_path.read_text())
    final_path = base/(name+'_epoch3.safetensors')
    assert sha(final_path) == final_receipt['checkpoint_sha256']
    final = load_file(str(final_path))
    assert set(final) == set(state['parameters'])
    for key,value in final.items():
        assert torch.equal(value,state['parameters'][key]),key
    optimizer_steps = sorted({int(row['step']) for row in state['optimizer']['state'].values()})
    assert optimizer_steps == [progress['updates']]
    result = {'passed':True,'time_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'fit':name,'completed_epochs':3,'committed_backprop_updates':progress['updates'],
              'execution_revision_sha256':revision_hash,
              'original_identity_fields_match_initial_metadata':True,
              'recovery_parameters_match_final_saved_model_exactly':True,
              'parameter_scalars':sum(value.numel() for value in final.values()),
              'all_Adam_parameter_steps_match_record':True,'schedule_digest_matches_final_record':True,
              'artifact_sha256':{p.name:sha(p) for p in [record_path,recovery,final_receipt_path,final_path]},
              'scope':'Completed recovery payload and final model/state identity. Saved effect/pair metric and neural-forward audits are separate; this does not establish predictive usefulness.',
              'code_sha256':sha(Path(__file__))}
    target = OUT/(name+'_completed_recovery_audit.json')
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['passed','fit','committed_backprop_updates','parameter_scalars']}),flush=True)


if __name__ == '__main__':
    main()

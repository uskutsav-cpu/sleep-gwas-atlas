"""Audit-only full-batch score-input resource probe with explicit provenance.

Training delegates to the actual matched trainer. Observation hooks retain its
forward arguments, preserve the initial atomic checkpoint, and run one complete
development-validation pass after the last audit checkpoint. No measured errors
are computed and no audit weights initialize scientific fits.
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
import time

os.environ.setdefault('PYTORCH_MPS_HIGH_WATERMARK_RATIO','0.4')
os.environ.setdefault('PYTORCH_MPS_LOW_WATERMARK_RATIO','0.3')
from permitted_AVI_experiment import OUT,ROOT,get_AVI_plan,sha


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--arm',choices=['scoreless','scores'],required=True)
    parser.add_argument('--audit-tag',required=True)
    args = parser.parse_args()
    assert args.audit_tag.startswith('AVI_resources_') and args.audit_tag.replace('_','').isalnum()
    processes = subprocess.run(['ps','-axo','pid,command'],check=True,capture_output=True,text=True).stdout
    forbidden = {'run_resumable_comparison.py','readonly_native_runtime.py','readonly_native_runtime_768.py',
                 'train_resumable_development.py','train_permitted_AVI_development.py','run_permitted_AVI_comparison.py'}
    for row in processes.splitlines()[1:]:
        fields = row.strip().split(None,1)
        if len(fields) != 2 or int(fields[0]) == os.getpid():
            continue
        parts = shlex.split(fields[1])
        if parts and Path(parts[0]).name.startswith('python'):
            assert not any(Path(p).name in forbidden for p in parts[1:]), 'No concurrent scientific fitting during resource audit.'
    plan = get_AVI_plan(for_execution_audit=True)
    original = importlib.import_module('train_permitted_AVI_development')
    torch,np = original.torch,original.np
    base_model = original.PermittedAVIModel

    class ObservedModel(base_model):
        def forward(self,*inputs,**kwargs):
            self._engineering_forward_inputs = inputs
            return super().forward(*inputs,**kwargs)

    original.PermittedAVIModel = ObservedModel
    directory = OUT/'AVI_recovery_audit_runs'/args.audit_tag
    directory.mkdir(parents=True,exist_ok=True)
    composition = {'time_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   'engineering_only':True,'PID':os.getpid(),'arm':args.arm,
                   'runtime_class':'ObservedModel inheriting the unchanged PermittedAVIModel forward',
                   'hooks':['Retain forward arguments without copying arrays.','Observe actual atomic recovery saves.','Full development inference after the committed update-1000 audit state.'],
                   'original_trainer_source_inventory_requires_this_addendum':True,
                   'weights_or_predictions_used_for_scientific_initialization_or_selection':False,
                   'code_sha256':{p.name:sha(p) for p in [Path(__file__),Path(original.__file__),
                                    Path(__file__).with_name('permitted_AVI_experiment.py'),
                                    Path(__file__).with_name('permitted_AVI_network.py'),
                                    Path(__file__).with_name('resumption.py')]}}
    composition_path = directory/'effective_runtime_composition.json'
    assert not composition_path.exists()
    composition_path.write_text(json.dumps(composition,indent=2)+'\n')
    real_save = original.save_recovery
    records = []
    started = time.monotonic()
    full_validation = None

    def state_digest(model):
        h = original.hashlib.sha256()
        for name,value in sorted(model.state_dict().items()):
            h.update(name.encode());h.update(value.detach().cpu().numpy().tobytes())
        return h.hexdigest()

    def observe(path,model,optimizer,device,identity,progress):
        nonlocal full_validation
        real_save(path,model,optimizer,device,identity,progress)
        if progress['updates'] == 0:
            snapshot = directory/'initial_recovery.pt'
            assert not snapshot.exists()
            os.link(path,snapshot)  # Future saves atomically replace the live path.
        if progress['updates'] == 1000:
            context = model._engineering_forward_inputs[1:]
            d,packed = context[:2]
            te = np.flatnonzero(d['split'] == 1)
            lengths = packed.length[d['bank_index'][te]]
            before = state_digest(model)
            cpu_rng,mps_rng = torch.get_rng_state(),torch.mps.get_rng_state()
            mode = model.training;model.eval();checked=0;last_print=time.monotonic()
            with torch.no_grad():
                for length in np.unique(lengths):
                    ids = te[lengths == length]
                    for offset in range(0,len(ids),256):
                        ix = ids[offset:offset+256]
                        predictions = model(ix,*context)
                        assert all(torch.isfinite(value).all() for value in predictions.values())
                        assert torch.all(predictions['magnitude'] >= 0)
                        checked += len(ix)
                        if time.monotonic()-last_print > 30:
                            assert shutil.disk_usage(ROOT).free > 300*1024**2, 'Stop the inference probe before exhausting disk reserve.'
                            print('AVI RESOURCE VALIDATION ROWS',checked,len(te),flush=True)
                            last_print=time.monotonic()
            model.train(mode)
            assert checked == len(te) and before == state_digest(model)
            assert torch.equal(cpu_rng,torch.get_rng_state()) and torch.equal(mps_rng,torch.mps.get_rng_state())
            full_validation = {'rows':checked,'native_lengths':np.unique(lengths).tolist(),
                               'prediction_indices_sha256':original.hashlib.sha256(te.tobytes()).hexdigest(),
                               'all_predictions_finite':True,'parameters_and_buffers_unchanged':True,
                               'CPU_and_MPS_RNG_unchanged':True,'experimental_errors_computed':False}
        records.append({'updates':progress['updates'],'phase':progress['phase'],
                        'elapsed_seconds':time.monotonic()-started,
                        'maximum_RSS_bytes_on_macOS':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                        'MPS_driver_allocated_bytes':torch.mps.driver_allocated_memory(),
                        'MPS_current_allocated_bytes':torch.mps.current_allocated_memory(),
                        'free_system_disk_bytes':shutil.disk_usage(ROOT).free})
        (directory/'resource_profile.json').write_text(json.dumps(records,indent=2)+'\n')
        assert records[-1]['free_system_disk_bytes'] > 300*1024**2, 'Stop with committed audit state before exhausting disk.'

    original.save_recovery = observe
    sys.argv = [str(original.__file__),'--recipe',plan['recipes'][0],'--held-study','none','--arm',args.arm,
                '--audit-tag',args.audit_tag,'--audit-stop-after-updates','1000']
    original.main()
    assert full_validation and records[-1]['updates'] == 1000
    result = {'passed':True,'actual_training_updates':1000,'arm':args.arm,
              'full_development_validation':full_validation,
              'maximum_RSS_bytes_on_macOS':max(row['maximum_RSS_bytes_on_macOS'] for row in records),
              'minimum_free_disk_at_recorded_checkpoints':min(row['free_system_disk_bytes'] for row in records),
              'resource_profile_sha256':sha(directory/'resource_profile.json'),
              'runtime_composition_sha256':sha(composition_path),
              'initial_recovery_sha256':sha(directory/'initial_recovery.pt'),
              'audit_weights_used_for_scientific_fits':False,'predictive_improvement_established':False}
    (directory/'resource_probe_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__ == '__main__':
    main()

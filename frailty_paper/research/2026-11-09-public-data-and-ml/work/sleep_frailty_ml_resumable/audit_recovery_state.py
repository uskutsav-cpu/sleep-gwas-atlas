"""Exercise the real fitter in separate processes: 4 steps versus 2 + reload + 2."""
import os
os.environ.setdefault('PYTORCH_MPS_HIGH_WATERMARK_RATIO','0.4')
os.environ.setdefault('PYTORCH_MPS_LOW_WATERMARK_RATIO','0.3')
import subprocess,time,shutil
from resumable_common import *
import torch


def run(tag,stop,resume=False):
    command=[sys.executable,'-u',str(Path(__file__).with_name('train_resumable_development.py')),
       '--recipe','adapt_last_stage','--held-study','none','--audit-tag',tag,'--audit-stop-after-updates',str(stop)]
    if resume:command.append('--resume')
    with (DEST/(tag+('_resumed' if resume else '')+'.log')).open('w') as handle:
        subprocess.run(command,stdout=handle,stderr=subprocess.STDOUT,check=True,
          env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'2','VECLIB_MAXIMUM_THREADS':'2'})
    print('COMPLETED RECOVERY AUDIT PROCESS',tag,stop,resume,flush=True)


def compare(a,b,path=''):
    if isinstance(a,torch.Tensor):
        assert isinstance(b,torch.Tensor) and a.dtype==b.dtype and torch.equal(a,b),path
        return a.numel()
    if isinstance(a,dict):
        assert set(a)==set(b),path
        return sum(compare(v,b[k],path+'/'+str(k)) for k,v in a.items())
    if isinstance(a,(list,tuple)):
        assert len(a)==len(b),path
        return sum(compare(v,z,path+'/'+str(i)) for i,(v,z) in enumerate(zip(a,b)))
    assert a==b,(path,a,b)
    return 0


def float32_replay(a,b,path,report):
    if isinstance(a,torch.Tensor):
        assert isinstance(b,torch.Tensor) and a.dtype==b.dtype and a.shape==b.shape,path
        if not a.is_floating_point():assert torch.equal(a,b),path
        else:
            assert torch.isfinite(a).all() and torch.isfinite(b).all()
            up=(torch.nextafter(a,torch.full_like(a,float('inf')))-a).abs()
            down=(torch.nextafter(a,torch.full_like(a,-float('inf')))-a).abs()
            unit=torch.maximum(up,down);gap=(a-b).abs()
            # A fixed eight representable float32 increments over four Adam
            # updates; direct checkpoint-state restoration remains exact.
            assert bool((gap<=8*unit).all()),(path,float(gap.max()))
            if not torch.equal(a,b):report.append({'tensor':path,'maximum_absolute_difference':float(gap.max()),
                'maximum_float32_ULP':float((gap/unit.clamp_min(torch.finfo(a.dtype).smallest_normal)).max())})
        return a.numel()
    if isinstance(a,dict):
        assert set(a)==set(b),path
        return sum(float32_replay(v,b[k],path+'/'+str(k),report) for k,v in a.items())
    if isinstance(a,(list,tuple)):
        assert len(a)==len(b),path
        return sum(float32_replay(v,z,path+'/'+str(i),report) for i,(v,z) in enumerate(zip(a,b)))
    assert a==b,(path,a,b)
    return 0


def main():
    start=time.monotonic();assert_versions_stable();get_stable_plan()
    assert not (DEST/'recovery_state_audit.json').exists()
    assert not (DEST/'recovery_audit_runs').exists(),'Inspect any partial audit before another attempt.'
    run('reference_four_steps',4);run('repeat_reference_four_steps',4);run('continued_four_steps',2)
    root=DEST/'recovery_audit_runs';suffix='development/adapt_last_stage_none_42_recovery.pt'
    partial=root/'continued_four_steps'/suffix
    shutil.copyfile(partial,DEST/'preserved_two_update_recovery_state.pt')
    partial.with_name(partial.name+'.pending').write_bytes(b'Uncommitted interrupted write; must not be loaded.')
    run('continued_four_steps',4,True)
    a=torch.load(root/'reference_four_steps'/suffix,map_location='cpu',weights_only=True)
    b=torch.load(partial,map_location='cpu',weights_only=True)
    report=[];parameters=float32_replay(a['parameters'],b['parameters'],'parameters',report)
    optimizer=float32_replay(a['optimizer'],b['optimizer'],'optimizer',report)
    fresh=torch.load(root/'repeat_reference_four_steps'/suffix,map_location='cpu',weights_only=True)
    fresh_report=[];float32_replay(a['parameters'],fresh['parameters'],'fresh_reference_parameters',fresh_report)
    float32_replay(a['optimizer'],fresh['optimizer'],'fresh_reference_optimizer',fresh_report)
    compare(a['CPU_rng'],b['CPU_rng'],'CPU_rng');compare(a['MPS_rng'],b['MPS_rng'],'MPS_rng')
    pa,pb=json.loads(a['progress_json']),json.loads(b['progress_json'])
    for key in ['updates','epoch','phase','schedule_sha256','NumPy_rng_json','curve']:compare(pa[key],pb[key],key)
    for key in ['objectives','anchors']:
        float32_replay(torch.tensor(pa[key],dtype=torch.float32),torch.tensor(pb[key],dtype=torch.float32),key,report)
    history=json.loads((root/'continued_four_steps/development/adapt_last_stage_none_42_resume_history.json').read_text())
    assert history[-1]['model_optimizer_and_CPU_MPS_random_state_restored_bit_for_bit']
    assert pa['updates']==pb['updates']==4
    write_json(DEST/'recovery_state_audit.json',{'passed':True,'actual_real_fitter_processes':4,'device':'mps',
       'uninterrupted_updates':4,'interrupted_prefix_updates':2,'resumed_total_updates':4,
       'parameter_tensor_scalars_compared':parameters,'optimizer_tensor_scalars_compared':optimizer,
       'initial_restored_model_optimizer_and_CPU_MPS_RNG_state_exact':True,
       'post_update_numeric_scope':'Final arithmetic checked within eight float32 representable increments over four updates. Bitwise identity of the full continued trajectory is not claimed. Initial restoration is checked exactly.',
       'post_update_differences':report,'fresh_uninterrupted_repeat_differences':fresh_report,
       'CPU_and_MPS_RNG_exactly_equal':True,'NumPy_sampler_schedule_exactly_equal':True,
       'uncommitted_pending_file_ignored':True,'engineering_models_not_used_to_initialize_scientific_fits':True,
       'new_study_outcomes_read':False,'predictive_improvement_established':False,'elapsed_seconds':time.monotonic()-start,
       'code_sha256':{p.name:sha(p) for p in Path(__file__).parent.glob('*.py')}})
    print('REAL DATA MPS RECOVERY STATE AUDIT PASS',flush=True)


if __name__=='__main__':main()

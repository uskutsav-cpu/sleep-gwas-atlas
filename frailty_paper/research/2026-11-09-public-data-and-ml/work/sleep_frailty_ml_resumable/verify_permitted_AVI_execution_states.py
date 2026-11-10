"""Verify actual AVI engineering states after all fitting workers have exited.

State comparisons are moved unchanged from the original coordinator so that its
Torch import does not remain resident beside each training/resource subprocess.
This process neither trains nor selects a model.
"""
import datetime,json,os,subprocess,shlex
from pathlib import Path
from permitted_AVI_experiment import OUT,AVI,PLAN_PATH,sha

def main():
    raw=subprocess.run(['ps','-axo','pid,command'],check=True,capture_output=True,text=True).stdout
    forbidden={'readonly_native_runtime.py','readonly_native_runtime_768.py','train_resumable_development.py','train_permitted_AVI_development.py','bench_permitted_AVI_resources.py'}
    for line in raw.splitlines()[1:]:
        fields=line.strip().split(None,1)
        if len(fields)!=2:continue
        args=shlex.split(fields[1])
        if args and Path(args[0]).name.startswith('python'):
            assert not any(Path(x).name in forbidden for x in args[1:]),'No training worker may overlap CPU state verification.'
    import torch
    from resumption import assert_exact
    from audit_readonly_runtime_recovery import numeric
    plan=json.loads(PLAN_PATH.read_text());recipe=plan['recipes'][0];runs=OUT/'AVI_recovery_audit_runs'
    def checkpoint(tag,arm):
        return runs/tag/'development'/f'{recipe}_none_{arm}_{plan["seed"]}_recovery.pt'
    def load(path):
        return torch.load(path,map_location='cpu',weights_only=True)

    continuations = {}
    direct_states = {}
    started_metadata = {}
    for arm in plan['arms']:
        direct_tag = f'AVI_execution_{arm}_direct4'
        continued_tag = f'AVI_execution_{arm}_continued4'
        prefix = checkpoint(continued_tag,arm)
        preserved = runs/continued_tag/'committed_two_update_prefix.pt'
        assert preserved.exists()  # Preserved by the trainer coordinator before continuation.
        assert json.loads(load(preserved)['progress_json'])['updates'] == 2
        a,b = load(checkpoint(direct_tag,arm)),load(checkpoint(continued_tag,arm))
        rows = []
        counts = {part:numeric(a[part],b[part],part,rows,1e-9,1e-6) for part in ['parameters','optimizer']}
        for part in ['CPU_rng','MPS_rng']:
            assert_exact(a[part],b[part],part)
        pa,pb = json.loads(a['progress_json']),json.loads(b['progress_json'])
        for part in ['updates','epoch','phase','schedule_sha256','NumPy_rng_json','curve']:
            assert_exact(pa[part],pb[part],part)
        for part in ['objectives','anchors']:
            numeric(torch.tensor(pa[part],dtype=torch.float32),torch.tensor(pb[part],dtype=torch.float32),part,rows,1e-9,1e-6)
        assert pa['updates'] == pb['updates'] == 4
        history_path = prefix.with_name(prefix.name.replace('_recovery.pt','_resume_history.json'))
        history = json.loads(history_path.read_text())
        assert history[-1]['updates'] == 2 and history[-1]['model_optimizer_and_CPU_MPS_random_state_restored_bit_for_bit']
        continuations[arm] = {'passed':True,'tensor_scalars':counts,'numeric_differences':rows,
                              'initial_restoration_exact':True,'RNG_and_schedule_exact':True,
                              'direct_sha256':sha(checkpoint(direct_tag,arm)),
                              'continued_sha256':sha(checkpoint(continued_tag,arm)),
                              'two_update_prefix_sha256':sha(preserved),'resume_history_sha256':sha(history_path)}
        direct_states[arm] = a
        started_path = checkpoint(direct_tag,arm).with_name(checkpoint(direct_tag,arm).name.replace('_recovery.pt','_started.json'))
        started_metadata[arm] = json.loads(started_path.read_text())
    common = ['recipe','seed','epochs','training_index_sha256','prediction_index_sha256','initial_parameter_sha256',
              'normalization','cell_medians','line_medians','line_parent','supported_cells','supported_lines',
              'context_transformations','AVI_scaling','AVI_bank_sha256','readonly_input_receipt','plan_sha256','code_sha256',
              'total_known_signed_training_rows','new_original_Siraj_signed_training_rows','new_original_neuronal_signed_training_rows']
    for key in common:
        assert started_metadata['scoreless'][key] == started_metadata['scores'][key],key
    for part in ['CPU_rng','MPS_rng']:
        assert_exact(direct_states['scoreless'][part],direct_states['scores'][part],part)
    progress = [json.loads(direct_states[arm]['progress_json']) for arm in plan['arms']]
    for part in ['updates','epoch','phase','schedule_sha256','NumPy_rng_json']:
        assert_exact(progress[0][part],progress[1][part],part)
    control = direct_states['scoreless']['parameters']['head.hidden.weight'][:,236:238]
    scored = direct_states['scores']['parameters']['head.hidden.weight'][:,236:238]
    assert torch.count_nonzero(control) == 0 and torch.count_nonzero(scored) > 0
    initials = [load(runs/f'AVI_resources_{arm}_v1'/'initial_recovery.pt') for arm in plan['arms']]
    for part in ['parameters','optimizer','CPU_rng','MPS_rng']:
        assert_exact(initials[0][part],initials[1][part],part)
    for arm in plan['arms']:
        assert started_metadata[arm]['plan_sha256']==sha(PLAN_PATH)
    result={'passed':True,'time_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'continuation_results':continuations,'all_common_initial_metadata_matched':True,
            'cross_arm_RNG_and_batch_schedule_exact':True,'score_columns_reached_training_gradients':True,
            'resource_probe_initial_states_equal_exactly':True,
            'comparison_code_moved_without_changing_numeric_tolerances':True,
            'fitting_worker_active_during_verification':False,'external_outcomes_read':False,
            'code_sha256':{p.name:sha(p) for p in [Path(__file__),Path(__file__).with_name('resumption.py'),Path(__file__).with_name('audit_readonly_runtime_recovery.py')]}}
    target=AVI/'matched_AVI_execution_state_verification.json';assert not target.exists()
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'passed':True,'verification_file_sha256':sha(target)}),flush=True)

if __name__=='__main__':main()

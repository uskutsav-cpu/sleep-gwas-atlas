"""Verify the actual attested wrapper before any scientific continuation."""
import datetime
import json
from pathlib import Path
import itertools

from audit_readonly_native_training import compare_states, RUNS, SUFFIX, torch
from audit_recovery_state import compare
from readonly_native_runtime import execution_revision
from readonly_native_inputs import sha
from resumable_common import DEST, ROOT


def numeric(a,b,path,rows,absolute,relative):
    if isinstance(a,torch.Tensor):
        assert isinstance(b,torch.Tensor) and a.shape == b.shape and a.dtype == b.dtype,path
        if a.is_floating_point():
            assert torch.isfinite(a).all() and torch.isfinite(b).all(),path
            gap = (a-b).abs()
            assert bool((gap <= absolute+relative*a.abs()).all()),(path,float(gap.max()))
            if not torch.equal(a,b):
                rows.append({'tensor':path,'scalars':a.numel(),'maximum_absolute_difference':float(gap.max())})
        else:
            assert torch.equal(a,b),path
        return a.numel()
    if isinstance(a,dict):
        assert set(a) == set(b),path
        return sum(numeric(v,b[k],path+'/'+str(k),rows,absolute,relative) for k,v in a.items())
    if isinstance(a,(list,tuple)):
        assert len(a) == len(b),path
        return sum(numeric(v,w,path+'/'+str(i),rows,absolute,relative) for i,(v,w) in enumerate(zip(a,b)))
    assert a == b,(path,a,b)
    return 0


def main():
    revision,revision_hash = execution_revision()
    direct = RUNS/'readonly_runtime_adapt_direct4'/SUFFIX
    continued = RUNS/'readonly_runtime_adapt_continued4'/SUFFIX
    reference = RUNS/'reference_four_steps'/SUFFIX
    protocol_path = DEST/'readonly_runtime_roundoff_control_protocol.json'
    protocol = json.loads(protocol_path.read_text())
    tolerance = protocol['additional_float32_comparison']
    absolute,relative = tolerance['absolute_tolerance'],tolerance['relative_tolerance']
    assert absolute == 1e-9 and relative == 1e-6
    names = ['reference_four_steps','repeat_reference_four_steps','readonly_native_adapt_direct4',
             'readonly_native_adapt_continued4','readonly_runtime_adapt_direct4',
             'readonly_runtime_adapt_continued4',*protocol['fixed_new_controls']]
    paths = {name:RUNS/name/SUFFIX for name in names}
    states = {name:torch.load(path,map_location='cpu',weights_only=True) for name,path in paths.items()}
    comparisons,legacy = {},{}
    for left,right in itertools.combinations(names,2):
        a,b = states[left],states[right]
        key = left+'_vs_'+right
        rows = []
        counts = {part:numeric(a[part],b[part],part,rows,absolute,relative) for part in ['parameters','optimizer']}
        for part in ['CPU_rng','MPS_rng']:
            compare(a[part],b[part],part)
        pa,pb = json.loads(a['progress_json']),json.loads(b['progress_json'])
        for part in ['updates','epoch','phase','schedule_sha256','NumPy_rng_json','curve']:
            compare(pa[part],pb[part],part)
        for part in ['objectives','anchors']:
            numeric(torch.tensor(pa[part],dtype=torch.float32),torch.tensor(pb[part],dtype=torch.float32),part,rows,absolute,relative)
        assert pa['updates'] == pb['updates'] == 4
        comparisons[key] = {'scalar_counts':counts,'differences':rows,'CPU_MPS_NumPy_RNG_and_schedule_exact':True}
        try:
            compare_states(paths[left],paths[right])
            legacy[key] = {'eight_ULP_passed':True}
        except AssertionError as error:
            legacy[key] = {'eight_ULP_passed':False,'first_failing_check':str(error)}
    # Keep the original same-wrapper continuation check; the new tolerance is
    # not used to hide a failed interruption/recovery arithmetic comparison.
    same_wrapper_legacy = compare_states(direct,continued)
    baseline = torch.load(reference,map_location='cpu',weights_only=True)
    original_identity = json.loads(baseline['identity_json'])
    for path in [direct,continued]:
        payload = torch.load(path,map_location='cpu',weights_only=True)
        identity = json.loads(payload['identity_json'])
        assert identity.pop('execution_revision_sha256') == revision_hash
        assert identity == original_identity
        addendum = path.with_name('adapt_last_stage_none_42_effective_runtime.json')
        assert json.loads(addendum.read_text())['execution_revision_sha256'] == revision_hash
    history_path = continued.with_name('adapt_last_stage_none_42_resume_history.json')
    history = json.loads(history_path.read_text())
    assert history[-1]['updates'] == 2 and history[-1]['model_optimizer_and_CPU_MPS_random_state_restored_bit_for_bit']
    result = {'passed':True,'time_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'execution_revision_sha256':revision_hash,'comparisons':comparisons,
              'additional_uniform_float32_tolerance':{'absolute':absolute,'relative':relative},
              'all_pair_legacy_eight_ULP_results':legacy,
              'same_wrapper_continuation_within_legacy_eight_ULP':same_wrapper_legacy,
              'initial_wrapper_vs_archived_reference_eight_ULP_check_failed':True,
              'initial_strict_failure_diagnostic_sha256':sha(DEST/'readonly_runtime_initial_eight_ULP_failure_diagnostic.json'),
              'control_protocol_sha256':sha(protocol_path),
              'all_original_and_wrapper_recovery_states_sha256':{name:sha(path) for name,path in paths.items()},
              'full_original_identity_retained_with_additional_execution_binding':True,
              'initial_parameter_optimizer_and_CPU_MPS_RNG_restoration_exact':True,
              'recovery_history_sha256':sha(history_path),
              'code_sha256':{str(p.relative_to(ROOT)):sha(p) for p in
                             [Path(__file__),Path(__file__).with_name('audit_readonly_native_training.py')]},
              'scope':'Actual attested execution: four updates and two plus two resumed, plus exactly two fresh unchanged-trainer controls. Report the original eight-ULP failure; all 28 state pairs must pass the uniform float32 comparison fixed before the two new controls. This is not bitwise full-trajectory identity or predictive validation. Scientific checkpoint migration still needs execution verification.',
              'scientific_predictive_success_established':False,
              'engineering_weights_used_for_scientific_initialization':False}
    path = DEST/'readonly_runtime_recovery_audit.json'
    assert not path.exists()
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'passed':True,'execution_revision_sha256':revision_hash,
                      'all_pairs_compared':len(comparisons),
                      'legacy_eight_ULP_failed_pairs':sum(not row['eight_ULP_passed'] for row in legacy.values()),
                      'maximum_absolute_difference':max((row['maximum_absolute_difference'] for r in comparisons.values() for row in r['differences']),default=0)}),flush=True)


if __name__ == '__main__':
    main()

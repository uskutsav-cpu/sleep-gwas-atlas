"""Verify the actual attested wrapper before any scientific continuation."""
import datetime
import json
from pathlib import Path

from audit_readonly_native_training import compare_states, RUNS, SUFFIX, torch
from readonly_native_runtime import execution_revision
from readonly_native_inputs import sha
from resumable_common import DEST, ROOT


def main():
    revision,revision_hash = execution_revision()
    direct = RUNS/'readonly_runtime_adapt_direct4'/SUFFIX
    continued = RUNS/'readonly_runtime_adapt_continued4'/SUFFIX
    reference = RUNS/'reference_four_steps'/SUFFIX
    comparisons = {'wrapper_vs_original':compare_states(reference,direct),
                   'wrapper_resumed_vs_direct':compare_states(direct,continued)}
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
              'full_original_identity_retained_with_additional_execution_binding':True,
              'initial_parameter_optimizer_and_CPU_MPS_RNG_restoration_exact':True,
              'recovery_history_sha256':sha(history_path),
              'code_sha256':{str(p.relative_to(ROOT)):sha(p) for p in
                             [Path(__file__),Path(__file__).with_name('audit_readonly_native_training.py')]},
              'scope':'Actual attested execution: four updates and two plus two resumed. Sustained allocation proof is separate. Migration of the interrupted scientific state must still be checked on execution.',
              'scientific_predictive_success_established':False,
              'engineering_weights_used_for_scientific_initialization':False}
    path = DEST/'readonly_runtime_recovery_audit.json'
    assert not path.exists()
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'passed':True,'execution_revision_sha256':revision_hash,
                      'max_observed_ULP':max((row['maximum_float32_ULP'] for r in comparisons.values() for row in r['differences']),default=0)}),flush=True)


if __name__ == '__main__':
    main()

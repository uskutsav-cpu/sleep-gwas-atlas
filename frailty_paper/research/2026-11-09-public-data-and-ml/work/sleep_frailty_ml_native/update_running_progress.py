"""Report actual process evidence, never infer liveness from a progress file."""
import subprocess,datetime
from native_common import *


def main():
    assert_versions_native();directory=DEST/'development';fits=[]
    for path in sorted(directory.glob('*_started.json')):
        r=json.loads(path.read_text());stem=path.name.removesuffix('_started.json');record=directory/(stem+'_record.json')
        if record.exists():
            done=json.loads(record.read_text());fits.append({'fit':stem,'state':'complete','completed_epochs':done['completed_epochs'],'updates':done['actual_backprop_updates']});continue
        p=subprocess.run(['ps','-p',str(r['PID']),'-o','command='],capture_output=True,text=True)
        live=p.returncode==0 and 'train_native_development.py' in p.stdout and '--recipe '+r['recipe'] in p.stdout and '--held-study '+r['held_study'] in p.stdout
        progress_path=directory/(stem+'_live_progress.json');progress={}
        if progress_path.exists():
            try:progress=json.loads(progress_path.read_text())
            except json.JSONDecodeError:progress={'file_write_in_progress':True}
        fits.append({'fit':stem,'state':'confirmed_running' if live else 'requires_terminal_state_inspection',
            'actual_PID_checked':r['PID'],'ps_command_matches_fit':live,'last_reported_progress':progress})
    audit=json.loads((DEST/'native_runtime_and_gradient_audit.json').read_text());plan=get_plan()
    write_json(DEST/'progress_and_readiness.json',{
       'checked_at_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
       'previous_goal_turn_classification':'Progress: finite representation comparison, convex optimality audit and original neuronal-label recovery completed.',
       'completed_now':['Native full-fragment model and matched comparison implemented',
          'Both arms share all documented signed labels, measured magnitudes, cell flags, initialization and batch-generation rules',
          'Independent actual-data finite differences verified gradients in two last-stage convolution tensors',
          'Actual MPS updates changed only the permitted last reader stage; prefix and BN running statistics preserved',
          'Size invariance and signed antisymmetry under allele exchange verified before and after updates'],
       'runtime_audit_passed':audit['passed'],'actual_process_checks':fits,
       'fixed_phase_one_roles':plan['phase_one_roles'],'three_epochs_per_arm_and_role':True,
       'phase_one_comparison_complete':(DEST/'phase_one_neural_gate.json').exists(),
       'further_source_fits_require_prespecified_neural_gate':True,
       'new_independent_test_prediction_or_confidence_success_established':False,
       'genetic_score_validation_or_patient_risk_model_established':False,
       'AlphaGenome_training_outputs_added':False,'prior_permissible_AVI_rejections_preserved':True,
       'new_Ulirsch_test_used_for_model_selection':False,'original_results_and_frozen_models_preserved':True,
       'scientific_goal_complete':False})
    (DEST/'PLAIN_LANGUAGE_PROGRESS.txt').write_text('''The matched DNA-model adaptation comparison is running.

Both versions use the same measured effects and documented direction labels.
The control keeps its DNA reader fixed; the other version can update its last
reading stage. Both see complete original DNA fragments at their actual
lengths. The model treats effect size consistently when the DNA versions are
exchanged. No AlphaGenome API-generated training targets or features were added.

Native-reader calculations, gradients and actual parameter updates passed
checks on real training examples. These are implementation checks, not proof
of useful biological prediction. Initial audited updates changed 290,896 scalars
in the permitted reader stage while the frozen prefix stayed unchanged.

The first comparison covers held regions and three held-out neuronal studies.
Each matched version has three fixed training passes. A predefined effect-size
and nearby/cell-difference gate decides whether to run the four remaining
held-study comparisons. Epochs and settings will not be extended to find a
favorable result. Earlier cropped-DNA core-tuning attempts are preserved.

Read progress_and_readiness.json for a fresh process check. A log or checkpoint
alone does not establish that a job is still running. Do not restart a fit
because observing its process timed out.

The reserved public blood-cell study remains out of model selection. Useful
effect sizes, cell differences, reliable confidence and a successful independent
added-input test remain required. No patient sleep/frailty prediction or
independent genetic-score validation is established.
''')
    print(json.dumps({'fits':fits,'scientific_goal_complete':False}),flush=True)


if __name__=='__main__':main()

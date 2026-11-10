"""Record actual live process evidence and scientific gaps without promotion."""
import subprocess,datetime
from stable_common import *


def main():
    assert_versions_stable();now=datetime.datetime.now(datetime.timezone.utc).isoformat()
    snapshot=subprocess.run(['ps','-axo','pid,etime,rss,command'],capture_output=True,text=True,check=True).stdout
    live=[r.strip() for r in snapshot.splitlines() if 'python' in r and any(s in r for s in ['audit_stable_runtime.py','train_stable_development.py','run_stable_comparison.py']) and '-c ' not in r]
    preflight=DEST/'stable_runtime_preflight.json';result=json.loads(preflight.read_text()) if preflight.exists() else None
    prior=json.loads((DEST/'benchmark_exposure/author_benchmark_exposure_receipt.json').read_text())
    final=json.loads((DEST/'benchmark_exposure/saturation_benchmark_component_receipt.json').read_text())
    records=list((DEST/'development').glob('*_record.json')) if (DEST/'development').exists() else []
    scientific={'effect_size_prediction':{'completed':False,'evidence':'Earlier completed model versions did not establish useful independent-study accuracy; this matched adaptation comparison remains unfinished.'},
       'nearby_and_cell_specific_effect_differences':{'completed':False,'evidence':'Existing fine cell-difference predictions remain weak; retain measured pair gates in the matched comparison.'},
       'useful_confidence_across_experiments':{'completed':False,'evidence':'Broad coverage or all-withheld outputs do not establish useful confidence. No newly promoted candidate.'},
       'AlphaGenome_added_input_untouched_test':{'completed':False,'evidence':'Earlier permitted AVI comparisons failed their gates. No ordinary restricted API outputs used for training. New model and permitted-input comparison must be fixed before revealing reserved outcomes.'}}
    report={'time_UTC':now,'authoritatively_live_process_rows':live,'CPU_runtime_preflight':result,
       'stable_scientific_fits_completed':len(records),'original_failed_run_preserved':True,
       'MPS_preflight_terminal_exit_code':1,'MPS_preflight_failure':'2.12 GiB allocated at a 2.13 GiB hard limit; no termination signal sent',
       'CPU_execution_plan_sha256':sha(DEST/'finite_stable_native_plan.json'),
       'new_test_exposure_screen':{'six_original_author_benchmark_files_verified':len(prior['benchmark_sources']),
            'eight_original_saturation_files_verified':len(final['published_saturation_benchmark_loci']),
            'candidate_designs_before':prior['candidate_designs_before'],'candidate_designs_after':final['candidate_designs_after'],
            'candidate_variants_after':final['candidate_variants_after'],'candidate_components_after':final['candidate_components_after'],
            'new_study_outcomes_used_for_selection':False,'remaining_supplementary_methods_review':True},
       'requirements':scientific,'candidate_promoted':False,'project_variants_authorized_for_ranking':False,'scientific_goal_complete':False,
       'next':'Poll the same confirmed-live CPU preflight session 23633 / PID 69451. Do not restart merely because it is slow. Require passed final preflight, inspect its actual evidence, then launch run_stable_comparison.py. Audit saved fits with audit_fitted_stable.py. Preserve original v1 and MPS failures.'}
    write_json(DEST/'current_progress_and_readiness.json',report)
    text=('The graphics-processor memory preflight failed at its hard allocation limit.\n'
          'A fresh CPU preflight uses the same measured effects and scientific model.\n'
          'Live processes verified at '+now+':\n'+'\n'.join(live)+'\n\n'
          'The scientific comparison has '+str(len(records))+' completed fits in this revision.\n'
          'Its prediction improvement is not established.\n\n'
          'We also screened the reserved study against six full original-reader benchmark\n'
          'files and eight original saturation experiments across four published loci.\n'
          'The conservative whole-component screen leaves 391 candidate fragment tests,\n'
          '134 DNA variants and 87 components, down from 1317 tests and 261 components.\n'
          'Supplementary-methods review remains before calling this a fully cleared test.\n'
          'No reserved outcome comparison or model selection used these measurements.\n\n'
          'Actual effect-size accuracy, fine cell differences, useful confidence and a\n'
          'successful permitted AlphaGenome added-input test all remain incomplete.\n')
    (DEST/'PLAIN_LANGUAGE_PROGRESS.txt').write_text(text)
    print(json.dumps({'live_processes':live,'preflight_complete':result is not None,'scientific_fits_completed':len(records),'scientific_goal_complete':False}),flush=True)


if __name__=='__main__':main()

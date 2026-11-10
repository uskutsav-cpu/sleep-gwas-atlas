"""Keep live execution evidence separate from prediction success."""
import datetime,subprocess
from bounded_common import *


def main():
    assert_versions_stable();now=datetime.datetime.now(datetime.timezone.utc).isoformat()
    snapshot=subprocess.run(['ps','-axo','pid,etime,rss,command'],capture_output=True,text=True,check=True).stdout
    live=[r.strip() for r in snapshot.splitlines() if 'python' in r and 'sleep_frailty_ml_bounded/' in r and any(s in r for s in ['audit_stable_runtime.py','train_stable_development.py','run_stable_comparison.py']) and '-c ' not in r]
    proof=DEST/'stable_runtime_preflight.json';partial=DEST/'stable_runtime_preflight_live.json'
    complete=json.loads(proof.read_text()) if proof.exists() else None
    latest=json.loads(partial.read_text())['observations'][-1] if partial.exists() else None
    equivalence=json.loads((DEST/'bounded_calculation_equivalence.json').read_text())
    avi=json.loads((DEST/'AVI_exposure/AVI_development_chromosome_exposure_receipt.json').read_text())
    development=DEST/'development';fits=list(development.glob('*_record.json')) if development.exists() else []
    assessment=DEST/'native_comparison_assessment.json'
    result={'time_UTC':now,'authoritatively_live_process_rows':live,'latest_runtime_observation':latest,
       'runtime_preflight_complete':complete is not None,'runtime_preflight_passed':complete['passed'] if complete else False,
       'calculation_and_training_gradient_equivalence_passed':equivalence['passed'],
       'scientific_fits_completed':len(fits),'scientific_comparison_assessment':json.loads(assessment.read_text()) if assessment.exists() else None,
       'candidate_promoted':False,'scientific_goal_complete':False,
       'AVI_comparison_candidate_pool':{'tests':avi['candidate_tests_after'],'variants':avi['candidate_variants_after'],
           'components':avi['candidate_components_after'],'outside_AVI_development_chromosomes':True,'already_sealed_count_QC_still_to_apply':True},
       'reserved_outcomes_used_for_model_selection':False,
       'remaining_requirements':['Useful independently tested actual effect-size accuracy','Useful nearby and cell-specific differences',
           'Empirically useful confidence across experiments','A successful matched permitted-AVI added-input comparison; no restricted API-output training'],
       'next':('Follow the currently verified comparison driver/worker; audit saved fits with audit_fitted_stable.py. Preserve incomplete fits before any restart.' if complete and complete['passed'] and live else
           'Inspect complete passed preflight criteria, then launch run_stable_comparison.py from this version.' if complete and complete['passed'] else
           'Poll confirmed-live runtime session 94677 / PID 70780 if still live. Require complete passed preflight and inspect every criterion before scientific fitting. Do not revive terminated CPU session 23633.')}
    write_json(DEST/'current_progress_and_readiness.json',result)
    text=('Same measured-effect prediction task; smaller native-reader pieces preserve every original scientific observation and loss.\n'
          'Independent checks reproduce original predictions and training gradients to numerical precision.\n'
          'This is a calculation check, not evidence of useful prediction accuracy.\n\n'
          'Verified live process state at '+now+':\n'+'\n'.join(live)+'\n'
          'Latest engineering update: '+str(latest['step'] if latest else 'not recorded')+' of 1000; a full validation pass also remains until the preflight is complete.\n'
          'Completed scientific fits in this revision: '+str(len(fits))+'.\n\n'
          'The reserved study has 391 candidate fragment tests after original-reader benchmark exclusions.\n'
          'A stricter matched AVI comparison pool has 183 candidate tests, 64 DNA changes and 47 linked groups outside AVI training/validation chromosomes.\n'
          'These counts precede the already sealed measurement-QC filter. The outcome results remain unused for model selection.\n\n'
          'Effect-size accuracy, cell and nearby differences, useful confidence and the permitted added-input test remain incomplete.\n')
    (DEST/'PLAIN_LANGUAGE_PROGRESS.txt').write_text(text)
    print(json.dumps({k:result[k] for k in ['authoritatively_live_process_rows','latest_runtime_observation','runtime_preflight_passed','scientific_fits_completed','scientific_goal_complete']}),flush=True)


if __name__=='__main__':main()

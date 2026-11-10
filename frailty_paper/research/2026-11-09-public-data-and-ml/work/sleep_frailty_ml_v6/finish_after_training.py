"""Serial verification, confidence fitting and one frozen whole-study test."""
import argparse,subprocess,sys,datetime
from common import *
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--bootstrap-pid',type=int,default=44410);a=parser.parse_args()
    print('Waiting on verified expanded-training bootstrap',a.bootstrap_pid,flush=True)
    while True:
        r=subprocess.run(['ps','-p',str(a.bootstrap_pid),'-o','command='],capture_output=True,text=True)
        if r.returncode or 'sleep_frailty_ml_v6/bootstrap_expanded_training.py' not in r.stdout:break
        time.sleep(10)
    assert (OUT/'final_fit_registry.json').exists(),'Training/refits did not finish; no outcome opening permitted.'
    started=(OUT/'final_test_evaluation_started.json').exists()
    steps=[] if started else ['audit_fitted_models.py','calibrate_uncertainty.py','check_scorer.py','prepare_reserved_features.py','freeze_final_model.py']
    steps+=['evaluate_reserved.py','additional_reserved_checks.py','write_final_report.py','build_delivery.py','audit_completion.py']
    for step in steps:
        path=WORK/step;assert path.exists(),step
        write_json(OUT/'finishing_progress.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'step':step,'reserved_outcomes_already_observed':(OUT/'final_test_evaluation_started.json').exists(),'scope':'Goal incomplete until final requirement audit and deliverable inspection.'});print('FINISHING STEP',step,flush=True)
        with (OUT/(path.stem+'.log')).open('a') as log:subprocess.run([sys.executable,'-u',str(path)],check=True,stdout=log,stderr=subprocess.STDOUT)
    print('ALL FIVE-STEP DELIVERABLES BUILT AND AUDITED; ROOT MUST INSPECT COMPLETION EVIDENCE BEFORE MARKING GOAL COMPLETE',flush=True)
if __name__=='__main__':main()

"""Sequential completion of the predeclared experiment; never restart a live fetch."""
import sys,time,subprocess,json,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).parent;DEST=ROOT/'outputs/sleep_frailty_ml_next';FINAL=DEST/'validation'

def run(script,*arguments):
    print('VALIDATION STEP START',script,*arguments,flush=True)
    subprocess.run([sys.executable,'-u',str(WORK/script),*arguments],check=True,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'2','VECLIB_MAXIMUM_THREADS':'2'})
    print('VALIDATION STEP COMPLETE',script,*arguments,flush=True)

def main():
    while not (DEST/'coverage_AVI_input_receipt.json').exists():
        progress=json.loads((DEST/'coverage_AVI_query_progress.json').read_text());pid=progress['PID'];process=subprocess.run(['ps','-p',str(pid),'-o','command='],capture_output=True,text=True)
        if process.returncode or 'fetch_coverage_AVI.py' not in process.stdout:
            if (DEST/'coverage_AVI_input_receipt.json').exists():break
            raise RuntimeError('The previously confirmed input fetch is no longer live and has no completion receipt. Inspect its session; do not silently restart it.')
        print('WAITING ON CONFIRMED LIVE AVI FETCH',pid,progress['windows_done'],'/',progress['windows_total'],flush=True);time.sleep(30)
    assert json.loads((DEST/'coverage_AVI_input_receipt.json').read_text())['passed']
    if not (FINAL/'AVI_fixed_execution_receipt.json').exists():run('run_fixed_coverage.py','--arm','AVI')
    if not (FINAL/'fixed_fits_independent_audit.json').exists():run('audit_fixed_coverage.py')
    if not (FINAL/'AVI_confidence_assessment.json').exists():run('calibrate_coverage_confidence.py','--arm','AVI')
    if not (FINAL/'AVI_confidence_independent_audit.json').exists():run('audit_coverage_confidence.py','--arm','AVI')
    if not (FINAL/'before_reserved_outcomes_lock.json').exists():
        run('check_pretest_inference_and_metrics.py');run('freeze_reserved_predictions.py')
    if not (FINAL/'reserved_outcome_evaluation.json').exists():run('evaluate_reserved_outcomes.py')
    if not (FINAL/'reserved_results_independent_audit.json').exists():run('audit_reserved_results.py')
    print('DECLARED COVERAGE EXPERIMENT COMPLETE; READ ACTUAL RESULTS BEFORE ANY SCIENTIFIC CLAIM',flush=True)

if __name__=='__main__':main()

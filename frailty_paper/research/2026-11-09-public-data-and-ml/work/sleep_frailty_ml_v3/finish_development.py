"""Finish locked selection/refits after the verified development job exits.

Stops before any reserved-test inference so the final artifact lock can be
inspected. It does not ask for user approval or redefine the scientific goal.
"""
from pathlib import Path
import sys,subprocess,time,json,argparse,hashlib
ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent;OUT=ROOT/'outputs/sleep_frailty_ai_v3'
def main():
    p=argparse.ArgumentParser();p.add_argument('--wait-for-cv',type=int,required=True);a=p.parse_args();print('Waiting for verified development job',a.wait_for_cv,flush=True)
    while True:
        s=subprocess.run(['ps','-p',str(a.wait_for_cv),'-o','command='],capture_output=True,text=True)
        if s.returncode or 'run_development_robustness.py' not in s.stdout:break
        time.sleep(10)
    assert (OUT/'development_robustness_finished.json').is_file()
    for file in ['summarize_development.py','train_conditional_stackers.py','select_final_recipes.py','final_fit.py','verify_final_fits.py']:
        print('RUN',file,flush=True);subprocess.run([sys.executable,'-u',str(WORK/file)],cwd=ROOT,check=True)
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();registry=json.loads((OUT/'final_model_registry.json').read_text());verification=json.loads((OUT/'final_fit_verification.json').read_text());assert verification['registry_sha256']==digest(OUT/'final_model_registry.json') and verification['no_reserved_test_scores_evaluated']
    lock={'status':'Final research artifacts selected, refitted and checked; reserved evaluation still pending','registry_sha256':digest(OUT/'final_model_registry.json'),'final_training_plan_sha256':digest(OUT/'final_training_locked_plan.json'),'verification_sha256':digest(OUT/'final_fit_verification.json'),'reserved_exposure_audit_sha256':digest(OUT/'reserved_exposure_audit.json'),'selected_role_names':registry['role_model_names'],'test_evaluated':False,'remaining':['One gated reserved Siraj chromosome evaluation with uncertainty and controls','Clear model card, limitations, study-specific validation and project-variant research predictions','Final reproducible package and evidence audit']};(OUT/'final_artifact_lock.json').write_text(json.dumps(lock,indent=2)+'\n');print('FINAL ARTIFACTS LOCKED; RESERVED EVALUATION PENDING',flush=True)
if __name__=='__main__':main()

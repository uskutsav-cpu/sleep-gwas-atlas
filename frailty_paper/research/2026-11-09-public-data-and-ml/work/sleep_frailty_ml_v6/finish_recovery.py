"""Finish the availability-driven recovery without ever reopening model selection."""
import subprocess,sys,time,datetime
from common import *
def main():
    deadline=time.monotonic()+1800
    while not (OUT/'additional_study_validation_lock.json').exists():
        if time.monotonic()>deadline:raise RuntimeError('Additional prediction lock not ready; no evaluation performed.')
        time.sleep(3)
    analysis=json.loads((OUT/'additional_analysis_code_lock.json').read_text())
    for p,h in analysis['files'].items():assert sha(p)==h,p
    original=json.loads((OUT/'final_model_lock.json').read_text())
    for p,h in {**original['files'],**original['model_checkpoints']}.items():assert sha(p)==h,p
    phases=[('evaluate_additional_study.py','reserved_study_results.json'),('additional_checks_v2.py','reserved_additional_checks.json'),('write_final_report.py','model_card.json'),('build_delivery.py','local_delivery_receipt.json'),('audit_completion.py','completion_audit.json')]
    for script,result in phases:
        if (OUT/result).exists():print('Existing completed phase retained',script,flush=True);continue
        write_json(OUT/'recovery_finishing_progress.json',{'phase':script,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'original_model_lock_sha256':sha(OUT/'final_model_lock.json')})
        print('START',script,flush=True)
        with (OUT/(script.removesuffix('.py')+'.log')).open('w') as log:subprocess.run([sys.executable,'-u',str(WORK/script)],stdout=log,stderr=subprocess.STDOUT,check=True,env=r_environment())
        assert (OUT/result).exists();print('DONE',script,flush=True)
    write_json(OUT/'recovery_finishing_progress.json',{'phase':'complete','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'completion_audit_sha256':sha(OUT/'completion_audit.json')});print('ALL FIVE MODEL WORKFLOW STEPS AUDITED AND DELIVERED; BIOLOGICAL WEAKNESSES DISCLOSED',flush=True)
if __name__=='__main__':main()

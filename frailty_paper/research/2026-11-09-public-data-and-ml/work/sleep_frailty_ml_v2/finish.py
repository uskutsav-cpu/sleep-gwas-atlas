"""Await specifically identified jobs, then verify and evaluate without refitting."""
from pathlib import Path
import argparse,subprocess,sys,time,json,hashlib
ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent;OUT=ROOT/'outputs/sleep_frailty_ai_v2'
p=argparse.ArgumentParser();p.add_argument('--student-pid',type=int,required=True);p.add_argument('--external-avi-pid',type=int,required=True);p.add_argument('--external-features-pid',type=int);a=p.parse_args()
for pid,expected in [(a.student_pid,'train_student.py'),(a.external_avi_pid,'fetch_avi.py'),(a.external_features_pid,'prepare_external_features.py')]:
    if pid is None:continue
    print('Waiting for verified job',pid,expected,flush=True)
    while True:
        state=subprocess.run(['ps','-p',str(pid),'-o','command='],capture_output=True,text=True)
        if state.returncode or expected not in state.stdout:break
        time.sleep(10)
    print('Job terminal',pid,flush=True)
assert json.loads((OUT/'student_status.json').read_text())['status'] in ['trained','not_triggered']
assert (WORK/'external_features_local.npz').is_file() and (OUT/'external_avi_retrieval_audit.json').is_file()
feature_receipt=json.loads((OUT/'external_feature_receipt.json').read_text())
assert feature_receipt['features_sha256']==hashlib.sha256((WORK/'external_features_local.npz').read_bytes()).hexdigest()
for script in ['evaluate.py','verify_metrics.py','score_project.py','report.py']:
    print('RUN',script,flush=True);subprocess.run([sys.executable,str(WORK/script)],cwd=ROOT,check=True)
(OUT/'EXPERIMENT_FINISHED.json').write_text(json.dumps({'state':'v2 finite search and external evaluation completed; broader strength/validation audit still required','goal_complete':False},indent=2)+'\n')
print('V2 EXPERIMENT FINISHED; inspect results and unresolved scientific validation',flush=True)

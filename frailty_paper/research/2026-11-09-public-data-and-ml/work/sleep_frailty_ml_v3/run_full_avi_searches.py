"""Wait on the actual annotated data job, then train; never unlock blind test."""
from pathlib import Path
import sys,subprocess,time,json,argparse
ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent;OUT=ROOT/'outputs/sleep_frailty_ai_v3'
p=argparse.ArgumentParser();p.add_argument('--avi-pid',type=int,required=True);p.add_argument('--biological-pid',type=int,required=True);p.add_argument('--prior-contrast-pid',type=int);a=p.parse_args()
jobs=[(a.avi_pid,'sleep_frailty_ml_v3/fetch_avi.py'),(a.biological_pid,'run_biological_searches.py')]
if a.prior_contrast_pid:jobs.append((a.prior_contrast_pid,'train_contrast.py --teacher prior'))
for pid,expected in jobs:
    print('Waiting for verified source job',pid,expected,flush=True)
    while True:
        state=subprocess.run(['ps','-p',str(pid),'-o','command='],capture_output=True,text=True)
        if state.returncode or expected not in state.stdout:break
        time.sleep(10)
    print('Source job terminal',pid,flush=True)
assert (OUT/'biological_searches_finished.json').is_file()
subprocess.run([sys.executable,'-u',str(WORK/'train_full_avi_classical.py')],cwd=ROOT,check=True)
for representation in ['hyena_chromatin','motif_chromatin']:
    for avi_input,task in [(False,'activity_avi'),(True,'activity_effect')]:
        cmd=[sys.executable,'-u',str(WORK/'train_biological_heads.py'),'--representation',representation,'--teacher','full','--tasks',task]
        if avi_input:cmd.append('--avi-input')
        print('RUN',representation,'full teacher','AVI input',avi_input,flush=True);subprocess.run(cmd,cwd=ROOT,check=True)
(OUT/'all_scheduled_searches_finished.json').write_text(json.dumps({'status':'All scheduled model fits complete; verification, robust selection, potential native biological fine-tuning and blind test remain','blind_test_evaluated':False,'goal_complete':False},indent=2)+'\n');print('ALL SCHEDULED SEARCHES FINISHED; BLIND TEST REMAINS LOCKED',flush=True)

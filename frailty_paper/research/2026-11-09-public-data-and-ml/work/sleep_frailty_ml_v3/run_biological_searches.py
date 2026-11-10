"""Run source-complete searches sequentially to bound local GPU memory."""
from pathlib import Path
import sys,subprocess,time,argparse,json
ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent;OUT=ROOT/'outputs/sleep_frailty_ai_v3'
p=argparse.ArgumentParser();p.add_argument('--wait-for-motif-pid',type=int);a=p.parse_args()
if a.wait_for_motif_pid:
    print('Waiting for verified motif extraction PID',a.wait_for_motif_pid,flush=True)
    while True:
        state=subprocess.run(['ps','-p',str(a.wait_for_motif_pid),'-o','command='],capture_output=True,text=True)
        if state.returncode or 'prepare_motif_features.py' not in state.stdout:break
        time.sleep(10)
assert (OUT/'motif_feature_receipt.json').is_file() and (OUT/'reference_chromatin_feature_receipt.json').is_file()
for name in ['hyena_chromatin','motif','motif_chromatin']:
    print('RUN BIOLOGICAL SEARCH',name,flush=True)
    subprocess.run([sys.executable,'-u',str(WORK/'train_biological_heads.py'),'--representation',name],cwd=ROOT,check=True)
(OUT/'biological_searches_finished.json').write_text(json.dumps({'status':'Three finite validation searches complete; broad-AVI comparison and blind test still pending','test_evaluated':False,'goal_complete':False},indent=2)+'\n');print('BIOLOGICAL SEARCHES FINISHED',flush=True)

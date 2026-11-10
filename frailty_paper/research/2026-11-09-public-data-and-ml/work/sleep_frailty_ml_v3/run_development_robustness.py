"""Finite, serial, development-only robustness comparisons after screening.

Never evaluate reserved chromosomes1–3. Checkpoint reload is independent of fitting.
"""
from pathlib import Path
import sys,subprocess,time,json,argparse,hashlib
ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent;OUT=ROOT/'outputs/sleep_frailty_ai_v3'
p=argparse.ArgumentParser();p.add_argument('--wait-for-full-contrast',type=int,required=True);a=p.parse_args()
print('Waiting for verified full contrast job',a.wait_for_full_contrast,flush=True)
while True:
    state=subprocess.run(['ps','-p',str(a.wait_for_full_contrast),'-o','command='],capture_output=True,text=True)
    if state.returncode or 'train_contrast.py --teacher full' not in state.stdout:break
    time.sleep(10)
assert (OUT/'contrast_full_val4-5_selection.json').is_file() and (OUT/'all_scheduled_searches_finished.json').is_file()
assert (OUT/'development_fold_sequence_audit.json').is_file()
def run(file,args=()):
    cmd=[sys.executable,'-u',str(WORK/file),*args];print('RUN',cmd,flush=True);subprocess.run(cmd,cwd=ROOT,check=True)
# All completed default-fold weights must reproduce saved predictions first.
run('check_validation_reload.py')
folds=[list(range(4,10)),list(range(10,16)),list(range(16,23))]
rules=json.loads((OUT/'development_robustness_prespecification.json').read_text());assert rules['folds']==folds
input_families=[]
for rep in ['hyena_chromatin','motif_chromatin']:
    path=OUT/(rep+'_full_avi_input_search_results.json');records=json.loads(path.read_text());input_families.extend([r for r in records if r['family']=='biological_seed_ensemble'])
input_winner=max(input_families,key=lambda r:r['validation']['macro_AP_skill']);input_rep=input_winner['features']
stage=json.loads((OUT/'stage1_search_results.json').read_text());classical=json.loads((OUT/'full_avi_classical_results.json').read_text())
pick=lambda rr:max(rr,key=lambda r:r['validation']['macro_AP_skill'])
baselines=[next(r for r in classical if r['family']=='direct_avi'),pick([r for r in stage if r['features']=='dna_pattern' and r['family'] in ['study_logistic','study_hist_boosting']]),pick([r for r in classical if r['features']=='chromatin']),pick([r for r in classical if r['features']=='dna_chromatin_avi'])]
locked={**rules,'screening_selected_AVI_input_ensemble':input_winner,'screening_selected_AVI_input_representation':input_rep,'classical_records':baselines,'data_sequence_audit_sha256':hashlib.sha256((OUT/'development_fold_sequence_audit.json').read_bytes()).hexdigest(),'status':'Locked before any three-fold comparison fitting','blind_test_evaluated':False}
path=OUT/'development_robustness_locked_config.json'
if path.exists():assert json.loads(path.read_text())==locked,'Do not overwrite a different CV plan'
else:path.write_text(json.dumps(locked,indent=2)+'\n')
done=[]
for fold in folds:
    val=','.join(map(str,fold));suffix='-'.join(map(str,fold))
    jobs=[('train_fold_baselines.py',['--validation-chromosomes',val],'development_baselines_val'+suffix),('train_contrast.py',['--teacher','full','--validation-chromosomes',val,'--tasks','effect,effect_avi'],'contrast_full_val'+suffix),('train_biological_heads.py',['--representation','motif_chromatin','--teacher','full','--validation-chromosomes',val,'--tasks','activity_effect,activity_avi'],'motif_chromatin_full_val'+suffix),('train_biological_heads.py',['--representation',input_rep,'--teacher','full','--avi-input','--validation-chromosomes',val,'--tasks','activity_effect'],input_rep+'_full_avi_input_val'+suffix)]
    for file,args,prefix in jobs:
        marker=OUT/(prefix+'_selection.json')
        if marker.is_file():print('Reuse completed family',prefix,flush=True)
        else:run(file,args)
        assert marker.is_file();done.append(prefix)
(OUT/'development_robustness_finished.json').write_text(json.dumps({'completed_families':done,'model_selection_and_blind_evaluation_pending':True,'blind_test_evaluated':False},indent=2)+'\n');print('DEVELOPMENT ROBUSTNESS FITS COMPLETE; BLIND TEST STILL LOCKED',flush=True)

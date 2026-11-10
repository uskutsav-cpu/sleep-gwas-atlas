"""Execute independent matched fits with bounded local process concurrency."""
from v9_common import *
import subprocess,time,concurrent.futures

FINAL=ROOT/'outputs/sleep_frailty_ai_v9_validation'


def main():
    lock=json.loads((FINAL/'final_settings_lock.json').read_text());directory=FINAL/'fits';directory.mkdir(exist_ok=True)
    tasks=[(s,a,seed) for s in lock['held_studies'] for seed in lock['seeds'] for a in ['sequence','AVI']]
    def run(task):
        study,arm,seed=task;name=f'fixed_{study}_{arm}_{seed}';record=directory/(name+'_record.json')
        if record.exists():
            r=json.loads(record.read_text());assert r['settings_lock_sha256']==sha(FINAL/'final_settings_lock.json');assert sha(Path(r['checkpoint']))==r['checkpoint_sha256'];assert sha(Path(r['predictions']))==r['predictions_sha256'];return name
        cmd=[sys.executable,'-u',str(WORK/'refit_fixed_candidate.py'),'--held-study',study,'--seed',str(seed)]+(['--avi'] if arm=='AVI' else [])
        print('START FIXED V9',name,flush=True)
        with (directory/(name+'.log')).open('w') as handle:subprocess.run(cmd,stdout=handle,stderr=subprocess.STDOUT,check=True,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'2','VECLIB_MAXIMUM_THREADS':'2'})
        print('DONE FIXED V9',name,flush=True);return name
    start=time.monotonic()
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:names=list(pool.map(run,tasks))
    for study in lock['held_studies']:
        for seed in lock['seeds']:
            a,b=[json.loads((directory/f'fixed_{study}_{arm}_{seed}_record.json').read_text()) for arm in ['sequence','AVI']]
            for key in ['shared_initial_parameter_sha256','training_index_sha256','prediction_index_sha256','completed_epochs']:
                assert a[key]==b[key],(study,seed,key)
    write_json(FINAL/'fixed_refit_execution_receipt.json',{'completed_fits':names,'count':len(names),'all_matched_arm_shared_weights_rows_and_epochs_pass':True,'elapsed_seconds':time.monotonic()-start,'fresh_model_performance_read':False})
    print('ALL FIXED V9 COMPLETE',len(names),flush=True)


if __name__=='__main__':main()

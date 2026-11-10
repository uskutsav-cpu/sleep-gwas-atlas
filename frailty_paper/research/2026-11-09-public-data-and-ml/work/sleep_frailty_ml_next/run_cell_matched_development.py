"""One finite input ablation before any newly reserved source is scored."""
import sys,subprocess,concurrent.futures,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *

DEST=ROOT/'outputs/sleep_frailty_ml_next'


def main():
    assert (DEST/'cell_matched_feature_receipt.json').exists();planpath=DEST/'finite_cell_matched_development_plan.json'
    plan={'recipes':['three_readers','WTC11_only'],'held_studies':['none','Salomon','Kosicki2025','Siraj','IGVF240K','Koesterich2023','McAfee2023','Retallick2026'],'epochs':8,'seed':42,
          'comparison':'Identical small head, shared initial weights, rows, losses, seeds and batch schedule. WTC11-only zeros the two added-reader features/readouts. Both use full actual source DNA and physical-length metadata.',
          'selection':'One common recipe and epoch by mean whole-source-omission raw magnitude MAE improvement over training-only cell/culture median. Report all source results and paired differences; no per-source selected epochs represented as validation.',
          'primary_priority':'Measured effect magnitude regression; paired nearby/cell difference and measured activity auxiliaries use the same previously audited losses. No detection labels or new-study outcomes.',
          'source_data':'Exactly the previously revealed V9 primary data; new blood and unused Siraj library outcomes are not included.',
          'scope':'Development only. Source eligibility, confidence transfer, matched AVI and a newly frozen external evaluation still required.',
          'stopping_rule':'Complete this finite feature comparison before selecting a follow-up. A good training/development result cannot authorize new-source predictions. Any outside evaluation requires a fresh settings, confidence and prediction lock.'}
    if planpath.exists():assert json.loads(planpath.read_text())==plan
    else:write_json(planpath,plan)
    output=DEST/'development';output.mkdir(exist_ok=True);tasks=[(recipe,source) for recipe in plan['recipes'] for source in plan['held_studies']]
    def run(task):
        recipe,source=task;name=f'{recipe}_{source}_{plan["seed"]}'
        if (output/(name+'_record.json')).exists():return name
        print('START CELL-MATCHED FIT',name,flush=True)
        with (output/(name+'.log')).open('w') as handle:subprocess.run([sys.executable,'-u',str(Path(__file__).with_name('train_cell_matched_development.py')),'--recipe',recipe,'--held-study',source],stdout=handle,stderr=subprocess.STDOUT,check=True,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'2','VECLIB_MAXIMUM_THREADS':'2'})
        print('DONE CELL-MATCHED FIT',name,flush=True);return name
    started=time.monotonic()
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:names=list(pool.map(run,tasks))
    for source in plan['held_studies']:
        a,b=[json.loads((output/f'{r}_{source}_42_record.json').read_text()) for r in plan['recipes']]
        for key in ['training_index_sha256','prediction_index_sha256','shared_initial_head_sha256','completed_epochs']:assert a[key]==b[key],(source,key)
    write_json(DEST/'cell_matched_development_execution_receipt.json',{'fits':names,'count':len(names),'matched_head_rows_and_initial_weights_verified':True,'new_source_model_performance_read':False,'elapsed_seconds':time.monotonic()-started});print('FINITE CELL-MATCHED DEVELOPMENT COMPLETE',flush=True)


if __name__=='__main__':main()

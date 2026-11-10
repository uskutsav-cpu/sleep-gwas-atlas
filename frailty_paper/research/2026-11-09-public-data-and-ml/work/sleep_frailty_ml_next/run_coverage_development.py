"""One declared expansion, one architecture, eight fixed epochs and seven omissions."""
import sys,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
DEST=ROOT/'outputs/sleep_frailty_ml_next'

def main():
    assert json.loads((DEST/'coverage_feature_receipt.json').read_text())['passed'];assert json.loads((DEST/'coverage_independent_data_audit.json').read_text())['passed']
    plan={'architecture':'Existing CellMatchedEffectHead, unchanged 66 features and 8,461 trainable head parameters; frozen three author readers.',
          'held_studies':['none','Salomon','Kosicki2025','Siraj','IGVF240K','Koesterich2023','McAfee2023','Retallick2026'],'epochs':8,'seed':42,
          'selection':'One common epoch maximizing mean seven-whole-source-omission MAE skill over training-only cell/culture median; tie takes earlier epoch. Held loci descriptive. No per-source best epoch reported as independent validation.',
          'expansion':'Exactly locked OL27 eligible components; all original array and partition prefixes preserved. Whole Siraj omission removes all its libraries, including old OL42 and new OL27.',
          'objective':'Same original activity_loss: main magnitude MAE independent of SE, plus fixed signed/activity and nearby/cell paired auxiliaries. No detection targets.',
          'additional_source_or_architecture_search_before_reserved_test':False,
          'after_selection':'Fix candidate epoch, run matched sequence/AVI seeds and source-omitted fits, calibrate and lock confidence plus predictions, then assess reserved blood and OL33 outcomes.',
          'reserved_outcomes_used':False,'data_sha256':sha(DEST/'coverage_training_arrays.npz'),'feature_receipt_sha256':sha(DEST/'coverage_feature_receipt.json')}
    path=DEST/'coverage_development_plan.json'
    if path.exists():assert json.loads(path.read_text())==plan
    else:write_json(path,plan)
    output=DEST/'coverage_development';output.mkdir(exist_ok=True);started=time.monotonic();names=[]
    for held in plan['held_studies']:
        name=f'three_readers_{held}_42';names.append(name)
        if (output/(name+'_record.json')).exists():continue
        print('START COVERAGE FIT',name,flush=True)
        with (output/(name+'.log')).open('w') as handle:subprocess.run([sys.executable,'-u',str(Path(__file__).with_name('train_coverage_development.py')),'--held-study',held],stdout=handle,stderr=subprocess.STDOUT,check=True,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'2','VECLIB_MAXIMUM_THREADS':'2'})
        print('DONE COVERAGE FIT',name,flush=True)
    curves={s:json.loads((output/f'three_readers_{s}_42_curve.json').read_text()) for s in plan['held_studies'][1:]};choices=[]
    for i in range(plan['epochs']):
        per={s:c[i]['macro_MAE_skill_vs_line_baseline_supported_cells'] for s,c in curves.items()};choices.append({'common_epoch':i+1,'mean_source_line_median_error_skill':float(np.mean(list(per.values()))),'per_source':per})
    selected=max(choices,key=lambda r:r['mean_source_line_median_error_skill']);local=json.loads((output/'three_readers_none_42_curve.json').read_text())[selected['common_epoch']-1]
    write_json(DEST/'coverage_development_assessment.json',{'selected':selected,'all_common_epochs':choices,'held_locus_metrics_at_common_epoch':local,'fits':names,
          'plan_sha256':sha(path),'all_fits_complete':True,'no_reserved_performance_calculated':True,'scientific_objective_complete':False,'no_additional_source_or_architecture_search_before_new_test':True,'elapsed_seconds':time.monotonic()-started})
    print('COVERAGE DEVELOPMENT COMPLETE',json.dumps(selected),flush=True)

if __name__=='__main__':main()

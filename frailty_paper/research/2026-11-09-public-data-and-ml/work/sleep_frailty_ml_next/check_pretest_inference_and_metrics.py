"""Actual saved-inference replay and a hand-computable macro/bootstrap check."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
from coverage_inputs import load_inputs
from coverage_validation_common import load_ensemble,predict_inputs,DEST,FINAL
from evaluate_reserved_outcomes import primary_metrics,pair_metrics

def main():
    assert not (FINAL/'before_reserved_outcomes_lock.json').exists();r,models,p=load_ensemble('none','sequence');d,pairs,bank,x=load_inputs();ids=p['observation_index'][:200];b=d['bank_index'][ids]
    q=predict_inputs(models,r,bank[b],x[b],d['cell'][ids].astype('int64'),d['line'][ids].astype('int64'),np.ones(len(ids),dtype=bool));errors={k:float(np.max(abs(q[k]-p[k][:len(ids)]))) for k in q};assert max(errors.values())<1e-6
    y=np.array([0.,1.,2.,4.,8.,12.]);pred=np.array([.2,1.2,2.2,3.,6.,9.]);strata=np.array(['a','a','a','b','b','b']);baseline=np.repeat(1.,6);groups=np.array(['chr1:0','chr1:1','chr1:2','chr1:0','chr1:1','chr1:2'])
    m=primary_metrics(y,pred,baseline,baseline,groups,strata,200,20261108);manual=float((sum(abs(y[:3]-pred[:3]))/3+sum(abs(y[3:]-pred[3:]))/3)/2);assert abs(m['macro_MAE']-manual)<1e-12;assert abs(m['macro_within_cell_line_Spearman']-1.)<1e-12
    pairs_check=pair_metrics(y,pred,baseline,np.array([[0,3],[1,4],[2,5]]),np.ones(6,dtype=bool),groups,strata,np.full(6,.2),{'factor':1.,'absolute_error_floor':.1},200,20261108)
    assert abs(pairs_check['difference_MAE']-2.2)<1e-12 and abs(pairs_check['difference_Spearman']-1.)<1e-12 and pairs_check['ordering_accuracy_nonzero_measured_differences']==1.
    record={'passed':True,'saved_sequence_ensemble_replayed_from_raw_features':len(ids),'maximum_errors':errors,'equal_line_macro_MAE_independently_checked':manual,'paired_difference_MAE_independently_checked':2.2,
            'cluster_bootstrap_smokecheck_repetitions':200,'reserved_inputs_or_outcomes_used':False,'code_sha256':{n:sha(Path(__file__).parent/n) for n in ['coverage_validation_common.py','evaluate_reserved_outcomes.py','check_pretest_inference_and_metrics.py']}}
    write_json(FINAL/'inference_and_metric_pretest_audit.json',record);print(json.dumps(record),flush=True)

if __name__=='__main__':main()

"""Measured differences on untouched NN-calibration regions; descriptive coverage."""
from final_validation_common import *
from evaluate_fresh_study import local_metrics,rho


def main():
    with np.load(DATA/'training_arrays.npz') as z:d={k:z[k] for k in z.files}
    with np.load(DATA/'measured_pairs.npz') as z:pairs={k:z[k] for k in z.files}
    results={}
    for arm in ['sequence','AVI']:
        with np.load(FINAL/f'{arm}_main_calibration_predictions.npz') as z:p={k:z[k] for k in z.files}
        indices=p['observation_index'];assert np.all(d['split'][indices]==2);inverse=np.full(len(d['effect']),-1,dtype='int32');inverse[indices]=np.arange(len(indices));arms={}
        cal=json.loads((FINAL/f'{arm}_confidence_calibration.json').read_text())['paired']
        for label,q in pairs.items():
            q=q[np.all(inverse[q]>=0,axis=1)];a,b=q.T;ia,ib=inverse[a],inverse[b];dy=d['effect'][a]-d['effect'][b];dp=p['magnitude'][ia]-p['magnitude'][ib]
            cell_baseline=p['cell_baseline'][ia]-p['cell_baseline'][ib];line_baseline=p['baseline'][ia]-p['baseline'][ib]
            keys=['|'.join(sorted([str(d['cell_names'][d['cell'][x]]),str(d['cell_names'][d['cell'][y]])])) for x,y in q]
            factors=np.array([cal[label]['by_cell_pair'].get(k,{'factor':cal[label]['global_factor']})['factor'] for k in keys]);half=factors*(p['raw_error_scale'][ia]+p['raw_error_scale'][ib])
            error=float(abs(dy-dp).mean());cellerr=float(abs(dy-cell_baseline).mean());lineerr=float(abs(dy-line_baseline).mean());strata=d['study'][a].astype('int64')*10000+d['cell'][a]*100+d['cell'][b];rs=[]
            for s in np.unique(strata):
                ix=strata==s
                if ix.sum()>=30:
                    value=rho(dy[ix],dp[ix])
                    if value is not None:rs.append(value)
            nonzero=dy!=0;answerable=(abs(dp)>half)&p['pretest_policy_accepted'][ia]&p['pretest_policy_accepted'][ib]
            arms[label]={'pairs':len(q),'regions':len(np.unique(d['group'][a])),'difference_MAE':error,'zero_difference_MAE':float(abs(dy).mean()),'cell_median_difference_MAE':cellerr,'line_median_difference_MAE':lineerr,
                         'skill_vs_zero_difference':1-error/float(abs(dy).mean()),'skill_vs_cell_median_difference':1-error/cellerr,'skill_vs_line_median_difference':1-error/lineerr,
                         'difference_Spearman':rho(dy,dp),'macro_within_fixed_study_cell_contrast_Spearman':float(np.mean(rs)) if rs else None,'fixed_contrast_strata':len(rs),
                         'ordering_accuracy_nonzero_measured_differences':float(np.mean(np.sign(dy[nonzero])==np.sign(dp[nonzero]))),
                         'direct_pair_interval_coverage_descriptive':float(np.mean(abs(dy-dp)<=half)),'pretest_answerable_pairs':int(answerable.sum()),'authorized_pairs':0,
                         'scope':'NN settings/fitting excluded these regions. Outcomes also informed source-omission interval scaling, so coverage here is descriptive, not independent uncertainty validation.'}
            np.savez_compressed(FINAL/f'{arm}_{label}_calibration_pair_diagnostics.npz',observation_pair=q,measured_difference=dy,predicted_difference=dp,cell_baseline_difference=cell_baseline,line_baseline_difference=line_baseline,interval_halfwidth=half,pretest_answerable=answerable)
        results[arm]=arms
    write_json(FINAL/'calibration_pair_comparisons.json',{'arms':results,'independent_new_cross_cell_replication':False,'project_or_cross_cell_authorized_answers':False})
    print('CALIBRATION PAIR COMPARISONS COMPLETE',flush=True)


if __name__=='__main__':main()

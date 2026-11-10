"""Check development measurement compatibility; never read reserved outcomes."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'sleep_frailty_ml_measurement'))
from measurement_common import *
from measurement_inputs import load_data
from scipy.stats import spearmanr

OUT=ROOT/'outputs/sleep_frailty_ml_direction'


def agreement(y, delta):
    good=np.isfinite(y)&np.isfinite(delta);y=y[good].astype('float64');delta=delta[good].astype('float64')
    if not len(y): return {'n':0}
    absolute=abs(delta);nonzero=y!=0
    return {'n':len(y),'magnitude_Spearman':float(spearmanr(y,absolute).statistic) if np.ptp(y)>0 and np.ptp(absolute)>0 else None,
       'median_observed_magnitude':float(np.median(y)),'median_absolute_activity_difference':float(np.median(absolute)),
       'mean_absolute_disagreement':float(abs(y-absolute).mean()),
       'median_absolute_disagreement':float(np.median(abs(y-absolute))),
       'fraction_absolute_disagreement_under_0_01':float(np.mean(abs(y-absolute)<.01)),
       'fraction_observed_at_least_twice_absolute_activity_difference':float(np.mean(y>=2*absolute))}


def main():
    assert_prior_versions();d,pairs,bank,x,extra=load_data(); result={}; train=d['split']==0
    active=np.isfinite(d['activity_reference'])&np.isfinite(d['activity_alternate'])
    delta=d['activity_alternate']-d['activity_reference']
    for code,name in enumerate(d['experiment_names']):
        allrows=(d['experiment']==code)&train
        if not allrows.any(): continue
        mask=allrows&active
        row={'training_rows':int(allrows.sum()),'with_two_measured_allele_activities':int(mask.sum()),
             'original_signed_targets_available':int((allrows&np.isfinite(d['signed_effect'])).sum()),
             'original_effect_vs_absolute_activity_difference':agreement(d['effect'][mask],delta[mask])}
        signed=allrows&active&np.isfinite(d['signed_effect'])
        if signed.any():
            a=d['signed_effect'][signed];b=delta[signed]
            row['known_signed_vs_activity_difference']={'n':int(signed.sum()),'max_absolute_disagreement':float(abs(a-b).max()),
                'signed_Spearman':float(spearmanr(a,b).statistic)}
        if mask.any():
            raw=bank[d['bank_index'][mask]]; per={}
            for reader in range(3):
                # Native bank coordinates 0/1 are REF/ALT activity readouts.
                predicted=raw[:,reader,1]-raw[:,reader,0]
                per[str(reader)]={'signed_activity_difference_Spearman':float(spearmanr(delta[mask],predicted).statistic),
                    'original_magnitude_Spearman':float(spearmanr(d['effect'][mask],abs(predicted)).statistic)}
            row['unmodified_reader_diagnostics']=per
        result[str(name)]=row
    write_json(OUT/'development_target_compatibility.json',{
       'development_training_partition_only':True,'reserved_outcomes_read':False,
       'data_sha256':sha(OLD/'coverage_training_arrays.npz'),'experiments':result,
       'RNA_activity_difference_not_assumed_equal_to_original_paired_GLM_effect':True,
       'sign_convention_inference_from_correlation_alone':False,
       'no_targets_or_frozen_models_changed':True})
    assert_prior_versions();print('DEVELOPMENT TARGET COMPATIBILITY AUDITED',len(result),flush=True)


if __name__=='__main__': main()

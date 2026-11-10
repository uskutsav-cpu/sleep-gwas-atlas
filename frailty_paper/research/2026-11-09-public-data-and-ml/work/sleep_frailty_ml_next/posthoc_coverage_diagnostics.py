"""Diagnose the revealed failed experiment; never relabel this as new validation."""
import sys,itertools
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
from scipy.stats import spearmanr
DEST=ROOT/'outputs/sleep_frailty_ml_next';FINAL=DEST/'validation'

def main():
    assert json.loads((FINAL/'reserved_results_independent_audit.json').read_text())['passed'];rows=list(csv_rows(FINAL/'reserved_joined_measurements_and_predictions.csv.gz'));lines=sorted({r['cell_line'] for r in rows if r['source']=='Siraj2026' and r['library']=='OL33'});by_line={};paired={}
    for line in lines:
        rr=[r for r in rows if r['source']=='Siraj2026' and r['library']=='OL33' and r['cell_line']==line];y=np.array([float(r['absolute_original_Log2Skew']) for r in rr]);p=np.array([float(r['sequence_magnitude']) for r in rr]);base=np.array([float(r['sequence_baseline']) for r in rr])
        by_line[line]={'measurements':len(rr),'observed_mean':float(y.mean()),'observed_median':float(np.median(y)),'training_only_typical_effect':float(base[0]),'prediction_mean':float(p.mean()),'prediction_median':float(np.median(p)),
                       'observed_IQR':float(np.quantile(y,.75)-np.quantile(y,.25)),'prediction_IQR':float(np.quantile(p,.75)-np.quantile(p,.25)),
                       'scope':'Observed spread includes measurement noise; compressed prediction spread alone does not prove biological underfitting.'}
    for first,second in itertools.combinations(lines,2):
        mapping={r['original_variant_GRCh37']:r for r in rows if r['source']=='Siraj2026' and r['library']=='OL33' and r['cell_line']==second};aa=[r for r in rows if r['source']=='Siraj2026' and r['library']=='OL33' and r['cell_line']==first];bb=[mapping[r['original_variant_GRCh37']] for r in aa]
        dy=np.array([float(a['absolute_original_Log2Skew'])-float(b['absolute_original_Log2Skew']) for a,b in zip(aa,bb)]);dp=np.array([float(a['sequence_magnitude'])-float(b['sequence_magnitude']) for a,b in zip(aa,bb)]);db=np.array([float(a['sequence_baseline'])-float(b['sequence_baseline']) for a,b in zip(aa,bb)])
        paired[first+'|'+second]={'designs':len(aa),'observed_mean_difference':float(dy.mean()),'observed_median_difference':float(np.median(dy)),'predicted_mean_difference':float(dp.mean()),'training_typical_effect_difference':float(db[0]),
                                 'model_difference_MAE':float(abs(dy-dp).mean()),'zero_difference_MAE':float(abs(dy).mean()),'training_typical_difference_MAE':float(abs(dy-db).mean()),'within_fixed_pair_Spearman':float(spearmanr(dy,dp).statistic)}
    result={'posthoc_now_revealed_data_only':True,'not_an_additional_independent_validation':True,'by_cell_line':by_line,'fixed_cell_pair_diagnostics':paired,
            'verified_failure':'The fresh independent-publication error does not improve on a simple baseline; nearby order is close to chance, and aggregate cross-cell order is below chance. Added AVI failed its locked retention rules. Both arms withheld every case.',
            'next_research_questions':['Separate publication/protocol offsets from variant-dependent errors using revealed data only.','Test whether retaining the pretrained readers individual channels, rather than scalar summaries, improves cell-specific learning in a finite development comparison. This is an untested hypothesis, not an explanation established by these results.','Audit experiment repeat/noise limits before interpreting observed variance as learnable biological signal.','Require new untouched suitable experimental data for any changed model; do not rescue the current model by retuning against these tests.'],
            'existing_models_or_acceptance_rules_changed':False,'prior_outcome_result_sha256':sha(FINAL/'reserved_outcome_evaluation.json'),'code_sha256':sha(Path(__file__))}
    write_json(FINAL/'posthoc_failure_diagnostics.json',result);print(json.dumps(result),flush=True)

if __name__=='__main__':main()

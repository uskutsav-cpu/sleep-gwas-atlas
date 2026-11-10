"""Post-evaluation interpretation only; no changes to frozen models or test rules."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
DEST=ROOT/'outputs/sleep_frailty_ml_next';FINAL=DEST/'validation';MS=ROOT/'outputs/nhanes_manuscript_audit'

def pct(x):return f'{100*x:.2f}%'

def main():
    assert json.loads((FINAL/'reserved_results_independent_audit.json').read_text())['passed'];assert not (FINAL/'paired_interpretation_integration_receipt.json').exists();r=json.loads((FINAL/'reserved_outcome_evaluation.json').read_text());p=r['paired_difference_results']['sequence']['OL33_three_known_lines'];local,cross=p['local'],p['cross'];ci=local['ordering_accuracy_95CI'];avi=r['primary_accuracy']['AVI']['blood_K562'];gain_ci=avi['gain_vs_sequence_95CI']
    heading='### Nearby changes and differences between cell lines'
    text=f'''{heading}

In the held library, the sequence model correctly ordered absolute measured effects in {pct(local['ordering_accuracy_nonzero_measured_differences'])} of {local['pairs']:,} nearby comparisons across {local['genomic_1Mb_windows']} one-megabase windows. The bootstrap interval for ordering was {pct(ci[0])} to {pct(ci[1])}, including chance ordering. Its paired effect-difference error improvement over predicting no nearby difference was {pct(local['relative_gain_vs_training_median_difference'])}, with an interval crossing zero. The independent blood publication supplied only two eligible nearby comparisons, which cannot establish local replication.

For {cross['pairs']:,} comparisons of the same physical designs between the three cell lines, ordering accuracy was {pct(cross['ordering_accuracy_nonzero_measured_differences'])}. The average correlation within fixed pairs of cell lines was {cross['macro_within_fixed_cell_pair_Spearman']:.4f}. Predicting no between-cell difference gave lower absolute error ({cross['zero_difference_MAE']:.6f}) than the model ({cross['difference_MAE']:.6f}). Improvement over the less accurate training-median difference therefore does not establish variant-specific cell learning. These results do not support cell-specific prioritization of sleep or frailty variants.

Adding AVI increased independent blood prediction error by {pct(-avi['relative_MAE_gain_vs_sequence'])} relative to sequence alone; the paired window bootstrap interval for the error increase was {pct(-gain_ci[1])} to {pct(-gain_ci[0])}. The small held-library gain did not satisfy the predefined retention rule. The composite input remains excluded from authorized biological predictions.

The nominal 90% intervals achieved high coverage with a median halfwidth of about 0.45 log2 units, while both arms withheld all cases. This is conservative withholding without useful answer availability. The accuracy requirements were a minimum 5% error improvement over both training-only typical-effect baselines, positive lower confidence bounds and a within-cell-line correlation of at least 0.1 with a positive lower bound. The added-input requirements also included at least 1% improvement in both primary tests and preservation of nearby ranking.
'''
    receipts=[]
    for path in [MS/'manuscript_with_nhanes.md',MS/'supplement/molecular_prediction_cell_matched.md']:
        before=path.read_text();assert heading not in before;old=sha(path);path.write_text(before+'\n'+text);assert path.read_text().startswith(before);receipts.append({'path':str(path),'entire_previous_prefix_preserved':True,'previous_sha256':old,'updated_sha256':sha(path),'previous_characters':len(before)})
    write_json(FINAL/'paired_interpretation_integration_receipt.json',{'files':receipts,'frozen_model_settings_or_predictions_changed':False,'verified_result_sha256':sha(FINAL/'reserved_outcome_evaluation.json'),'code_sha256':sha(Path(__file__))})
    (FINAL/'figures').mkdir(exist_ok=True)
    (FINAL/'figures/figure_legend.txt').write_text('Reserved laboratory DNA-effect tests. Error improvement is relative to the training-only cell/culture median; positive values mean lower error. Whiskers show 95% paired one-megabase-window bootstrap intervals. The dotted line marks the prespecified 5% minimum improvement. Coverage is shown against a nominal 90% target. Answer availability is the fraction accepted by the input policy before targets were joined; it is zero for both arms. The blood K562 experiment is an independent publication; OL33 is a held library within the existing Siraj publication/lab. No patient-level score validation is shown.\n')
    print(text,flush=True)

if __name__=='__main__':main()

"""Separate scientific findings from implementation completion and prior evidence."""
from final_validation_common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def percent(value):return f'{100*value:.2f}%'


def main():
    result=json.loads((FINAL/'fresh_outcome_evaluation.json').read_text());selection=json.loads((FINAL/'common_epoch_selection.json').read_text());paired=json.loads((FINAL/'calibration_pair_comparisons.json').read_text())
    seq,avi=[result['arms'][a] for a in ['sequence','AVI']];sm,am=[a['all_cases'] for a in [seq,avi]];comparison=result['AVI_comparison'];gate=seq['full_model_quality_gate_passed'];delivery=seq['uncertainty_and_delivery_gate_passed']
    rows=[]
    for arm in ['sequence','AVI']:
        a=result['arms'][arm];m=a['all_cases'];rows.append({'arm':arm,**{k:v for k,v in m.items() if k!='bootstrap'},'all_case_interval_coverage':a['all_case_interval_coverage'],
                         'pretest_policy_accepted_cases':a['pretest_policy_accepted_cases'],'authorized_cases':a['authorized_cases'],'full_model_quality_gate_passed':a['full_model_quality_gate_passed'],'uncertainty_and_delivery_gate_passed':a['uncertainty_and_delivery_gate_passed']})
    write_csv(FINAL/'fresh_study_summary.csv',rows)
    report=f'''SLEEP–FRAILTY MOLECULAR MODEL: V9 RESULTS

What was actually predicted
The size of a DNA change's effect in a laboratory reporter experiment, in original log2 units. This is not a person's risk of insomnia or frailty.

Development
613,161 measured effects. Earlier McAfee and Retallick outside results were explicitly reassigned to development; their original adverse V7 results are preserved. The selected frozen-reader model trained for five epochs. Its mean development improvement over a cell/culture-specific median was {percent(selection['selection']['mean_source_skill'])}; selection used seven whole-study omissions and is not independent replication.

New study, with predictions and confidence settings frozen first
610 eligible baseline effects from the journal's original Matoba et al. table. The source includes assay-active elements in primary neural progenitors from one male donor. Independence from the earlier primary-cell donor has not been established. Cases include detected and undetected allelic effects; no significance filtering was applied.

Sequence-only error: {sm['MAE']:.6f}.
Training-only cell/culture median error: {sm['line_or_cell_baseline_MAE']:.6f}.
Improvement over that median: {percent(sm['relative_MAE_improvement_vs_line_baseline'])}.
Effect-size rank correlation: {sm['magnitude_Spearman']:.4f}.
Full usefulness gate passed: {gate}.

Confidence
Nominal 90% interval coverage in the new study: {percent(seq['all_case_interval_coverage'])}.
Cases accepted by the input policy frozen before the test: {seq['pretest_policy_accepted_cases']} of 610.
Cases authorized after all outside-performance protections: {seq['authorized_cases']}.
Uncertainty and delivery gate passed: {delivery}.
Calibration used whole-study-omitted errors, with direct calibration of paired differences. It gives no universal cross-study guarantee. Broad intervals or withholding every answer are not useful predictive success.

Added AVI score
Error with AVI: {am['MAE']:.6f}.
Improvement over the matched sequence-only model: {percent(comparison['relative_MAE_improvement_vs_sequence'])}.
Retain AVI under all frozen gates: {comparison['retain_AVI']}.
Nearby comparisons: {comparison['nearby_pairs']} pairs in {comparison['nearby_1Mb_windows']} one-megabase windows. Sufficient local confirmation: {comparison['local_confirmation_sufficient']}.
AVI is a composite containing AlphaGenome-derived and other information. This cannot isolate a benefit of AlphaGenome itself. Original V7 rejection remains unchanged.

Cell differences and nearby variants
Held calibration-region results are in calibration_pair_comparisons.json. Their confidence coverage is descriptive because these outcomes also informed interval scaling through source-omitted predictions. The new experiment tests one cell context, so it cannot independently validate cross-cell interactions. Project native genomic contexts have no verified reporter activity and receive no authorized answers from this active-element-only validation.

Still missing
Independent human genetic-score validation; a same-participant sleep/fatigue-excluded frailty GWAS; stronger variant-to-gene evidence; independent relevant-cell validation of differences; source assay raw-QC reconstruction; full pretrained genomic-neighborhood exclusion. NHANES remains separate observational evidence and does not resolve those genetic gaps.

Changing the model after seeing this test would make this study development data for that changed model. Another untouched experiment would then be required.
'''
    (FINAL/'PLAIN_LANGUAGE_RESULTS.txt').write_text(report)
    selected_arm='AVI' if comparison['retain_AVI'] else 'sequence'
    useful=bool(result['arms'][selected_arm]['uncertainty_and_delivery_gate_passed'])
    audit={'requirements':{
        '1_predict_actual_effect_size':{'implementation':'completed','useful_independent_prediction_established':bool(gate),'evidence':'fresh_outcome_evaluation.json: arms.sequence.all_cases and full_model_quality_gate_passed'},
        '2_learn_cell_and_nearby_differences':{'implementation':'paired losses and direct paired-difference interval calibration completed','independent_cross_cell_validation_established':False,'fresh_nearby_evidence':seq['nearby_comparisons'],'evidence':'calibration_pair_comparisons.json and fresh_outcome_evaluation.json'},
        '3_fix_confidence_across_experiments':{'implementation':'whole-source-omitted error models, region-balanced calibration, unfamiliar-input withholding completed','useful_independent_confident_prediction_established':bool(delivery),'evidence':'sequence_confidence_calibration.json and fresh_outcome_evaluation.json'},
        '4_test_AlphaGenome_as_added_input':{'matched_composite_test_completed':True,'retain_composite':bool(comparison['retain_AVI']),'AlphaGenome_specific_benefit_established':False,'evidence':'fixed_refit_execution_receipt.json and fresh_outcome_evaluation.json: AVI_comparison'}},
        'scientific_objective_complete':False,'reason':'Training and validation are implemented, but useful cross-cell molecular prediction and independent relevant-cell validation remain unproven. Implementation completion is not the requested scientific success.',
        'Nature_Communications_genetic_manuscript_ready':False,'human_disease_prediction_or_genetic_score_validation':False,'prior_scientific_analyses_preserved':True}
    write_json(FINAL/'requirement_completion_audit.json',audit)
    write_json(FINAL/'model_card.json',{'intended_use':'Research diagnostics for absolute reporter-effect magnitude; no person-level disease risk prediction.',
               'training_observations':613161,'contributing_magnitude_studies':7,'unused_metadata_study_label':'Kircher2019 has zero primary magnitude observations.','final_fit_regions':'Original training plus validation roles; calibration and fresh-source overlap excluded.',
               'pretraining':'Real author WTC11 MPRA-LegNet weights. Fixed convolutional reader; trainable culture/cell activity and magnitude heads.',
               'new_external_study':result['scope'],'useful_confident_prediction_established':useful,'authorized_project_contexts':[],
               'AVI_retained':bool(comparison['retain_AVI']),'validation_limitations':result['not_established'],
               'release_warning':'Raw outputs are diagnostics. Authoritative authorization fields must remain false wherever the external gates fail; no raw score sorting presented as biological ranking.'})
    # Preserve the entire old manuscript verbatim and append the new molecular
    # result as an explicitly separate, later benchmark.
    manuscript=ROOT/'outputs/nhanes_manuscript_audit/manuscript_with_nhanes.md';old=manuscript.read_text();title='## Coverage-amended molecular prediction benchmark (V9)';assert title not in old
    snapshot=FINAL/'manuscript_before_V9.md';snapshot.write_text(old)
    appendix=f'''

{title}

A later, explicitly developmental amendment reused the already revealed McAfee and Retallick measurements, preserving the original V7 results above. The V9 magnitude dataset contained 613,161 measurements. One line-aware architecture, with either fixed or partially adapted pretrained representations, was compared using a single common epoch across seven whole-study omissions. The selected fixed-reader, five-epoch candidate improved development mean absolute error over a training-only cell/culture median by {percent(selection['selection']['mean_source_skill'])}; this selected development result is not independent evidence of accuracy. Three-seed final and source-omitted ensembles were trained with and without AVI using matched rows, epochs and shared initial weights. Independent numerical replay and fitted-gradient checks verified computation, not biological usefulness.

Before joining new outcomes, models, error predictors, empirical 90% interval calibration, acceptance rules and all new-study predictions were frozen by cryptographic hash. The new experiment comprised 610 eligible original absolute fitted log2 coefficients from [Matoba et al. (2026)](https://doi.org/10.1038/s41593-026-02454-2), after source, construct orientation, reciprocal genome-coordinate and overlap checks. The assay-active source subset includes detected and undetected allelic effects; there was no significance filter. Its repeated cultures derive from one primary neural-progenitor donor, whose independence from the earlier primary-cell donor is not established. Original barcode QC and coefficient fitting were not independently reconstructed, signed source-allele conventions remain unverified, and full pretrained genomic-neighborhood independence is not established.

Sequence-only absolute error was {sm['MAE']:.6f}, compared with {sm['line_or_cell_baseline_MAE']:.6f} for the training-only cell/culture median ({percent(sm['relative_MAE_improvement_vs_line_baseline'])} relative improvement; negative values mean worse prediction). Effect-magnitude Spearman correlation was {sm['magnitude_Spearman']:.4f}. The frozen full-model usefulness gate passed: **{gate}**. Nominal 90% intervals covered {percent(seq['all_case_interval_coverage'])} of new measurements; {seq['pretest_policy_accepted_cases']} cases passed the input policy before outcomes were joined, and {seq['authorized_cases']} remained authorized after outside-performance protection. These results do not establish a universal confidence guarantee.

AVI changed error by {percent(comparison['relative_MAE_improvement_vs_sequence'])} relative to the matched sequence-only model. All retention gates passed: **{comparison['retain_AVI']}**. Local confirmation used {comparison['nearby_pairs']} nearby pairs in {comparison['nearby_1Mb_windows']} one-megabase windows; sufficient local confirmation was **{comparison['local_confirmation_sufficient']}**. Because AVI combines AlphaGenome-derived and other predictors, this experiment cannot attribute a benefit uniquely to AlphaGenome. The earlier V7 retention decision remains unchanged.

Direct paired-difference results on held NN-calibration regions compare predictions with zero, cell-median and culture-median differences. Those outcomes also informed source-omission interval scaling, so their coverage is descriptive. The new one-cell experiment cannot independently establish variant-specific differences between cell types. No project genomic context receives an authorized answer from this assay-active-only validation. This molecular work does not validate insomnia/frailty scores in people, establish genetic correlation or causation, or strengthen target-gene evidence by itself. The NHANES sensitivity result and all original genetic analyses remain unchanged. Detailed methods, complete adverse comparisons and readiness gaps are in the V9 molecular supplement.
'''
    manuscript.write_text(old+appendix);assert manuscript.read_text().startswith(old)
    supplement=ROOT/'outputs/nhanes_manuscript_audit/supplement/molecular_prediction_v9.md'
    supplement.write_text('# Coverage-amended laboratory-effect prediction supplement\n\n'+report+'\n\nSource-linked tables and reproducible code are in `../../sleep_frailty_ai_v9_validation/` and `../../../work/sleep_frailty_ml_v9/`. The new result tables are `fresh_study_summary.csv`, `fresh_study_predictions_and_observations.csv`, and `calibration_pair_comparisons.json`. This is molecular evidence, separate from NHANES and the original genetic analyses.\n')
    write_json(FINAL/'manuscript_integration_receipt.json',{'original_prefix_sha256':hashlib.sha256(old.encode()).hexdigest(),'snapshot':str(snapshot),'new_manuscript_sha256':sha(manuscript),'original_entire_prefix_preserved':True,'new_supplement':str(supplement),'observational_and_genetic_results_modified':False})
    source_rows=list(csv_rows(FINAL/'fresh_study_predictions_and_observations.csv'));y=np.array([float(r['measured_absolute_original_log2_effect']) for r in source_rows]);fig,axes=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    dev=selection['selection']['per_source_skill'];ax=axes[0,0];ax.barh(list(dev),[100*v for v in dev.values()],color=['#3577a5' if v>=0 else '#b15943' for v in dev.values()]);ax.axvline(0,color='black',lw=.8);ax.set(xlabel='Error reduction over cell/culture median (%)',title='A  Development studies (selected settings)')
    ax=axes[0,1];labels=['Cell median','Culture median','Sequence','With AVI'];values=[sm['cell_baseline_MAE'],sm['line_or_cell_baseline_MAE'],sm['MAE'],am['MAE']];bars=ax.bar(labels,values,color=['#999999','#555555','#3577a5','#8c60a7']);ax.bar_label(bars,fmt='%.3f',padding=3);ax.set(ylabel='Mean absolute error (log2 units)',title='B  New study: 610 measured effects');ax.set_ylim(0,max(values)*1.2)
    ax=axes[1,0];p=np.array([float(r['sequence_magnitude']) for r in source_rows]);ax.scatter(y,p,s=9,alpha=.4,color='#3577a5');limit=max(y.max(),p.max());ax.plot([0,limit],[0,limit],color='#555555',ls='--',lw=1);ax.set(xlabel='Measured effect magnitude',ylabel='Predicted effect magnitude',title=f'C  New-study ranking: Spearman {sm["magnitude_Spearman"]:.3f}',xlim=(0,limit*1.03),ylim=(0,limit*1.03))
    ax=axes[1,1];labels=['90% interval\ncoverage','Input-policy\nacceptance','Authorized\nanswers'];values=[seq['all_case_interval_coverage'],seq['pretest_policy_accepted_cases']/610,seq['authorized_cases']/610];bars=ax.bar(labels,np.array(values)*100,color=['#3577a5','#cfac52','#555555']);ax.bar_label(bars,fmt='%.1f%%',padding=3);ax.axhline(90,color='#777777',ls='--',lw=.8);ax.set(ylabel='New-study measurements (%)',ylim=(0,110),title='D  Coverage and answer availability')
    path=FINAL/'figures';path.mkdir(exist_ok=True);fig.savefig(path/'development_and_fresh_study_results.png',dpi=200);plt.close(fig)
    assert_previous_sources();print(report,flush=True)


if __name__=='__main__':main()

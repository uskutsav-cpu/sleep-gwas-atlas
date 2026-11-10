"""Evidence-based V8 progress; never equate training completion with success."""
from datetime import datetime, timezone
from v8_common import *


def main():
    assert_previous_intact()
    data = json.loads((OUT/'measured_training_data_audit.json').read_text())
    enrichment = json.loads((OUT/'measured_activity_uncertainty_enrichment_audit.json').read_text())
    fresh = json.loads((OUT/'fresh_matoba_physical_design_audit.json').read_text())
    publisher = json.loads((OUT/'fresh_matoba_publisher_provenance_addendum.json').read_text())
    avi = json.loads((OUT/'AVI_comparison_input_receipt.json').read_text())
    records = {}
    for path in OUT.glob('*_record.json'):
        if path.name.startswith('._'):
            continue
        r = json.loads(path.read_text())
        assert r['actual_backprop_updates'] > 0 and sha(r['checkpoint']) == r['checkpoint_sha256']
        records[path.name.removesuffix('_record.json')] = r
    transfer = {}
    for recipe in ['frozen', 'pooled_tune']:
        curves = []
        for study in ['Kosicki2025', 'IGVF240K', 'Salomon', 'Siraj']:
            curves.append(json.loads((OUT/f'{recipe}_{study}_42_curve.json').read_text()))
        macro = np.mean([[r['macro_MAE_skill_supported_cells'] for r in c] for c in curves], axis=0)
        e = int(np.argmax(macro))
        transfer[recipe] = {'best_common_epoch': e+1, 'mean_development_whole_study_error_skill_at_that_epoch': float(macro[e]), 'all_epoch_values': macro.tolist(), 'not_independent_final_validation': True}
    status = {
        'generated_UTC': datetime.now(timezone.utc).isoformat(),
        'completed_measured_data_fits': len(records),
        'actual_backprop_updates_across_fits': sum(r['actual_backprop_updates'] for r in records.values()),
        'small_pretrained_model_actual_parameters': json.loads((OUT/'small_model_runtime_and_gradient_audit.json').read_text())['total_parameter_scalars'],
        'prepared_magnitude_measurements': data['observations'],
        'main_fit_training_measurements': records['pooled_tune_none_42']['training_rows'],
        'fresh_test_candidate_measured_effects': fresh['reserved_pairs_with_source_original_baseline_coefficient'],
        'journal_baseline_equivalence_verified': publisher['journal_baseline_measurement_equivalence_established'],
        'additional_measured_allele_activities': enrichment['reconstructed_Kosicki_allele_activity_rows'],
        'all_measured_allele_activity_pairs_now_available': enrichment['after_paired_activity_rows'],
        'measured_standard_error_rows_now_available': enrichment['empirical_three_replicate_paired_mean_SE_rows']+enrichment['original_author_GLM_SE_rows'],
        'common_epoch_whole_study_transfer': transfer,
        'withheld_region_results': {k: r['best']['metrics'] for k, r in records.items() if r['held_study'] == 'none'},
        'numerical_correctness_verified': json.loads((OUT/'fitted_effect_network_independent_audit.json').read_text())['passed'],
        'fresh_AVI_composite_inputs_prepared': avi['available'],
        'fresh_outside_study_model_performance_calculated': False,
        'AVI_matched_training_and_retention_decision_completed': False,
        'new_model_cross_study_confidence_validated': False,
        'strong_prediction_or_reliable_cell_specific_ranking_established': False,
        'human_genetic_score_validation_established': False,
        'entire_project_scientifically_finished': False,
        'original_NHANES_genetics_and_V3_V6_V7_evidence_preserved': True,
        'next_required_work': [
            'Use genuinely measured per-allele activities and measurement uncertainty in an explicitly documented development comparison; do not simply add epochs to the current weak transfer recipe.',
            'Judge actual variant-specific cell effects against both zero and cell-median differences, including within-cell-pair rankings.',
            'Choose a recipe using development sources only; fit uncertainty and a matched AVI-added arm with training-only transformations.',
            'Lock all settings, weights, test predictions and retention rules before calculating the fresh study performance.',
            'Report adverse results and usable-answer availability; preserve the broader genetic and disease-validation gaps.'
        ],
    }
    write_json(OUT/'V8_RESEARCH_PROGRESS.json', status)
    fine = records['end_to_end_none_42']['best_development_score']*100
    text = f'''V8 SLEEP–FRAILTY MODEL: RESEARCH PROGRESS

As of {status['generated_UTC']}.

We completed {len(records)} measured-data development fits using a small pretrained DNA-activity model with 857,825 parameters. Existing pretrained weights were genuinely adjusted, including the last convolutional stage in the larger fit. Independent NumPy arithmetic and fitted-gradient checks passed. These checks establish that the training calculation works, not that predictions are biologically useful.

The prepared dataset contains {data['observations']:,} effect-size measurements. The main fit used {status['main_fit_training_measurements']:,} training observations, with separate genomic regions for model selection and calibration. The larger fit improved average error by about {fine:.2f}% on held development regions. Study-to-study transfer is still weak: the pooled recipes do not consistently beat a typical-effect guess across the source-omission comparisons. These are development comparisons used for choosing settings, not independent final validation. Much of the apparent between-cell correlation is explained by average differences between cell types; variant-specific cell learning remains weak.

We recovered another 20,725 real paired-allele activity measurements and now have {enrichment['after_paired_activity_rows']:,} such pairs. We also recovered measured uncertainty for {status['measured_standard_error_rows_now_available']:,} effects, keeping original GLM errors distinct from empirical errors across three assay replicates. Missing errors stay missing. The original training arrays and completed fits were preserved.

A free public brain-cell study supplies 610 candidate final-test effects after physical-DNA, coordinate, duplication, pretraining-context and prior-experiment overlap checks. The journal's original baseline tables were checked against the archived data and agree numerically. Both detected and undetected allelic effects are retained. The data cover active DNA regions in cells from one male donor, not independent people, older adults or a representative genome-wide sample. No project candidate variant has a matching measured effect in this reserved dataset.

Public permissible AVI composite inputs are available for all 610 candidate test alleles. This prepares a future matched comparison; it is not training on experimental truth from AlphaGenome and does not reverse the previous V7 decision to exclude AVI.

The next development experiment should learn from both measured DNA versions and their uncertainty. Confidence checks, the matched AVI experiment and the locked new-study test remain unfinished. The model is not a validated sleep/frailty risk model, genetic-score validation or proof of causal genes. No scientific completion is claimed.
'''
    (OUT/'V8_RESEARCH_PROGRESS.txt').write_text(text)
    print('V8 EVIDENCE INVENTORY WRITTEN', len(records), status['actual_backprop_updates_across_fits'], flush=True)


if __name__ == '__main__':
    main()

"""Fix final molecular tests before observing any fresh-study model result."""
from v8_common import *


def main():
    path = OUT/'V8_external_validation_protocol.json'
    protocol = {
        'fresh_study': 'Matoba et al. 2026, original journal Supplementary Table 3b and physically verified GEO/author constructs',
        'fresh_subset': 'All 610 eligible original baseline fitted coefficients, including detected and undetected allelic effects. Availability and physical/source audits determine eligibility; no model or outcome-based hit selection.',
        'fresh_genome_design_sha256': sha(DATA/'reserved_matoba_physical_design.csv'),
        'official_journal_source_sha256': sha(SOURCES/'Matoba_journal_Supplementary_Table3.xlsx'),
        'magnitude_target': 'Absolute original fitted log2 allelic coefficient. Do not infer its signed convention from genomic reference labels. Do not substitute activity differences or AlphaGenome predictions as test truth.',
        'limitations': ['Original assay-active subset, not a representative all-variant sample', 'Primary neural progenitors from one male donor, repeated cultures, not independent people', 'Source raw barcode QC and coefficient estimation not independently reconstructed', 'Exact pretrained context matches excluded; full pretrained genomic-neighborhood exclusion not established', 'No project allele has a matching measured coefficient in this fresh subset', 'No human insomnia/frailty genetic-score validation or causal-gene test'],
        'development_and_locks': {
            'model_selection': 'Development source/locus comparisons only; one common recipe and epoch. All adverse comparisons retained.',
            'main_final_fit': 'Refit on original training+validation regions only, with all new-test overlapping whole components excluded; calibration regions remain outside fitting.',
            'ensemble_seeds': [42, 73, 101],
            'matched_AVI': 'Identical source rows, target values, epoch count, seeds, shared initial weights and batch curriculum. AVI is appended as raw score, phred and availability after training-only scaling; added input weights start at zero.',
            'before_fresh_performance': ['Freeze selected architecture, epochs and source data', 'Freeze all final and source-omitted fitted weights', 'Freeze uncertainty and abstention policy', 'Save both arms predictions, intervals and eligibility before joining source coefficients'],
            'future_model_changes': 'A changed model cannot reuse this revealed study as an untouched confirmatory test.'
        },
        'primary_metrics': {
            'all_case_magnitude_MAE': 'Mean absolute error in raw original log2 units; paired comparison with training-only cell median and between model arms.',
            'ranking': 'Effect-magnitude Spearman globally and within genomic regions; report nearby pairs separately.',
            'bootstrap': '2,000 paired resamples of 1-Mb genomic windows, seed 20261108; state that windows approximate dependence and are not ancestry-specific LD blocks.',
            'full_model_usefulness_gate': {'relative_MAE_improvement_over_training_cell_median_at_least': .05, 'paired_MAE_improvement_CI95_lower_above_zero': True, 'magnitude_Spearman_at_least': .10, 'ranking_CI95_lower_above_zero': True},
            'cell_differences': 'On the untouched calibration loci, compare actual matched-cell differences against zero and training-cell median differences; also evaluate within fixed study/cell contrast. These loci are independent of model selection; interval calibration on them is not an independent confidence-coverage test.',
            'nearby_ordering': 'Same source experiment and genomic component, within 25 bp. No significance or effect-gap screen. Report unanswered comparisons and uncertainty; individual interval overlap is not a calibrated paired-order guarantee.'
        },
        'uncertainty': {
            'target': .90,
            'error_predictor_training': 'Predictions made while a whole development study and its genomic components were omitted; fit error predictors using only original training-role observations of those predictions. Calibration-role observations do not choose neural settings.',
            'error_models': {'family': 'HistGradientBoostingRegressor', 'max_iter': 100, 'max_leaf_nodes': 7, 'min_samples_leaf': 100, 'l2_regularization': 10, 'early_stopping': False, 'random_state': 42, 'targets': ['actual absolute model error, squared-error mean regression', 'actual absolute cell-median baseline error, squared-error mean regression', '0.90 quantile of actual absolute model error']},
            'allowed_inputs': 'Predicted magnitude, signed contrast, mean allele activity, global/local representation summaries, cell/family and number of same-cell training studies. No measured effect, activity, SE, p value or significance label is an inference input.',
            'calibration': 'Fixed-source-omission predictions on original calibration regions; region-balanced residual/error-scale ratios. Report both observation-weighted and region-balanced coverage and interval width. No universal cross-study guarantee.',
            'abstention_input_policy': {'known_modeled_cell': True, 'same_cell_training_source_count_at_least': 1, 'maximum_absolute_standardized_error_feature': 10, 'predicted_relative_error_improvement_over_cell_median_at_least': .05, 'predicted_absolute_error_improvement_at_least_log2_units': .005, 'supported_conditions': ['baseline', 'untreated', 'unreported']},
            'external_delivery_protection': {'frozen_policy_accepted_cases_at_least': 30, 'accepted_genomic_windows_at_least': 20, 'all_case_interval_coverage_at_least': .90, 'accepted_case_interval_coverage_at_least': .90, 'accepted_error_improvement_CI95_lower_above_zero': True, 'positive_effect_ranking': True, 'all_case_full_model_usefulness_gate_must_also_pass': True},
            'not_success': ['Broad intervals alone', 'High coverage obtained by withholding every answer', 'A large number of code checks', 'Post hoc restriction to favorable variants'],
        },
        'AVI_retention_gate': {
            'relative_MAE_improvement_over_sequence_only_at_least': .01,
            'paired_genomic_window_MAE_improvement_CI95_lower_above_zero': True,
            'local_ranking_point_drop_at_most': .01,
            'report_local_group_count_and_low_power': True,
            'if_insufficient_local_groups': 'Report insufficient local confirmation; do not claim a safe nearby-ranking improvement. Retention does not establish useful predictions unless the separate quality gate passes.',
            'source': 'Permissive downloaded composite AVI only. Any gain cannot be attributed uniquely to AlphaGenome.',
            'original_V7_AVI_rejection_unchanged': True,
        },
        'scope': 'Predeclared research validation and conditional release protection. Does not promise successful prediction, general coverage, disease causation or a Nature Communications-ready genetic manuscript.',
    }
    if path.exists():
        assert json.loads(path.read_text()) == protocol
    else:
        write_json(path, protocol)
    print('EXTERNAL VALIDATION AND RETENTION RULES FIXED BEFORE FRESH PERFORMANCE', flush=True)


if __name__ == '__main__':
    main()

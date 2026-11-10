"""Report frozen adverse results clearly; no model, decision or threshold changes."""
from common import *


def main():
    immutable = [
        'final_model_and_prediction_lock.json', 'reserved_study_results.json',
        'additional_neural_results.json', 'additional_neural_prediction_lock.json',
        'combined_model_recommendation.json', 'inference_registry.json',
    ]
    before = {n: sha(OUT/n) for n in immutable}
    primary = json.loads((OUT/'reserved_study_results.json').read_text())
    extra = json.loads((OUT/'additional_neural_results.json').read_text())
    combined = json.loads((OUT/'combined_model_recommendation.json').read_text())
    zero = json.loads((OUT/'posthoc_zero_cell_difference_comparator.json').read_text())
    domain = json.loads((OUT/'posthoc_context_length_diagnostic.json').read_text())
    neuron = extra['conditions']['baseline']['sequence_only']
    assert combined['all_cell_groups_withheld'] and not combined['retain_AVI_added_input']
    assessment = {
        'useful_cross_study_effect_magnitude_prediction_established': False,
        'reliable_variant_specific_cell_differences_established': False,
        'authorized_cell_groups': [],
        'additional_neuron_MAE_relative_to_training_cell_median': neuron['MAE']/neuron['baseline_MAE']-1,
        'primary_AVI_relative_error_improvement': primary['AVI_comparison']['relative_macro_MAE_improvement'],
        'additional_neuron_AVI_relative_error_improvement': extra['AVI_baseline_confirmation']['relative_MAE_improvement'],
        'AVI_retained_after_both_checks': False,
        'AVI_interpretation': 'Absolute error improved in both outside studies. The second study failed the predeclared local-ranking tolerance; only four eligible local groups contributed, so this is not strong evidence of biological harm. AVI combines AlphaGenome and other predictors.',
        'uncertainty_interpretation': 'Observed 97–99% coverage accompanies broad ranges. This does not establish useful precision or a universal coverage guarantee. Suppressing all answers is not predictive success.',
        'posthoc_zero_cell_difference_comparator': zero,
        'posthoc_domain_diagnostic': {
            '150_vs_151bp_length_failure_supported': False,
            'primary_fraction_beyond_domain_limit': domain['datasets']['primary_native150']['fraction_above_domain_limit12'],
            'reason': 'Only 0.12% of primary physical inputs exceed the feature-domain limit. All primary rows fail the frozen multi-study cell-support requirement; this is not evidence of a universal input-length artifact.',
            'scope': 'Diagnostic after outside results; no retuning or new confirmatory test.'
        },
        'frozen_input_policy_vs_final_delivery': 'The original policy accepted all 547 baseline neuronal measurements and none of the primary-study measurements. Subsequent predeclared outside-performance protection withholds all cell groups; the original measurements and predictions remain unchanged.',
    }
    write_json(OUT/'scientific_failure_assessment.json', assessment)
    plan = {
        'status': 'Proposed next experiments, not completed or proven improvements',
        'priority_order': [
            {'priority': 1, 'change': 'Learn measured activity of both DNA versions and their signed difference before taking its magnitude.', 'why': 'The current magnitude model barely transfers and does not learn reliable variant-specific cell differences.', 'guardrail': 'Confirm original allele sign conventions and physical assay orientation; keep incompatible experimental units separate.'},
            {'priority': 2, 'change': 'Benchmark a small MPRA-trained network, then compare frozen weights with partial fine-tuning.', 'candidates': ['MPRA-LegNet', 'Malinois for contexts with compatible source assays', 'PARM only for promoter-compatible loci'], 'why': 'The current generic DNA encoder was fixed; only small prediction heads learned measured outcomes.', 'guardrail': 'Audit pretraining overlap and licences; activity accuracy is not variant-effect accuracy. Published full MPRA-LegNet training requires more RAM than this 8 GB Mac; do not promise local full-scale training.'},
            {'priority': 3, 'change': 'Add real matched experiments across relevant cell types, including unchanged or near-zero effects.', 'why': 'All primary test rows lacked support from two development studies for their cell group; a cell label alone cannot supply missing biology.', 'guardrail': 'Do not duplicate the same measurements or replace experimental targets with AlphaGenome predictions. Tumour cells and immature neurons do not establish older-adult disease biology.'},
            {'priority': 4, 'change': 'Require improvement over both the training-cell median and the no-cell-difference comparator on development studies left out in turn.', 'why': 'The apparent between-cell improvement disappears against a zero-difference comparator.', 'guardrail': 'Report genomic-region splits, full-study transfer, effect error, nearby ordering, usable answer fraction and interval width together.'},
            {'priority': 5, 'change': 'Reserve a fresh study before model changes and use it once after locking the new model.', 'candidate': 'Matoba et al., GSE289858', 'availability': 'Public source barcode-to-150-bp physical-DNA table reached successfully; source schema checked. Full genome, orientation, duplicate, donor and neighborhood audits remain.', 'guardrail': 'The two revealed V7 outside studies are benchmarks, not untouched tests for a revised model. No Matoba model performance has been calculated.'},
        ],
        'sources': [
            'https://github.com/autosome-ru/human_legnet',
            'https://github.com/vansteensellab/PARM',
            'https://pmc.ncbi.nlm.nih.gov/articles/PMC10441439/',
            'https://www.nature.com/articles/s41593-026-02454-2',
            'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289858',
        ],
        'disease_validation_separate': 'None of these molecular experiments establishes a person-level sleep/frailty prediction model, independent genetic-score replication, disease causation or the missing sleep/exhaustion-excluded frailty GWAS.'
    }
    write_json(OUT/'prioritized_model_improvement_plan.json', plan)
    report = '''SLEEP–FRAILTY MODEL: VERIFIED RESULTS AND NEXT STEPS

The trained model is not strong enough to rank the project's DNA changes reliably. In the separate neuron study, its average error was about 7.2% higher than simply guessing the typical effect for that cell group. The other outside study also showed weak prediction. No cell group passed all checks for releasing rankings.

What we trained
Small prediction networks learned from measured laboratory outcomes using backpropagation. The underlying pretrained DNA reader stayed fixed. Development contained 605,885 measurements of actual effect magnitude plus 27,191 auxiliary measurements kept in a separate unit. These are measurements, not independent people or unique DNA changes.

What the outside tests showed
The primary test contained 9,602 measurements across two cell types from 4,820 verified DNA pairs. The additional baseline neuron test contained 547 measurements. Neither study was used to choose or fit the prediction networks. Original source estimates were checked, but source raw barcode processing was not independently reconstructed.

Adding AVI, a score combining AlphaGenome with other information, reduced average error by 2.71% in the primary study and 1.06% in the neuron study. It failed the second study's previously fixed check on ordering nearby DNA changes. Only four nearby groups contributed to that check, so we do not claim strong evidence of biological harm. The combined rule therefore excludes AVI. Its improvement cannot be assigned to AlphaGenome alone.

Predictions of differences between cells were also weak. A later descriptive check found that assuming no difference performed better than either model. We preserved the original results and labeled this extra comparison as added after seeing the tests.

Wide prediction ranges covered 97–99% of observed effects. Broad ranges and withholding every answer do not make this a useful model. The original input rules accepted all 547 neuron measurements; the final outside-performance check still withheld this cell group because prediction was worse than the simple comparison.

Project coverage
The supplement includes 498 distinct variants: the 488 molecular candidates combined with the score markers. It preserves 1,351 gene-evidence rows. Single-letter changes have 470 supported inputs; 28 insertions/deletions remain unsupported. All released answer and rank fields are blank. Separate raw projections are unvalidated diagnostics, not evidence of harmlessness, disease mechanisms or causal genes.

How to improve it
First learn the measured activity of both DNA versions and their difference. Compare a small model already trained on these laboratory experiments with the current fixed DNA reader. Add real examples tested in several relevant cell types, and require it to beat simple guesses when an entire study is left out. A new public brain-cell study, GSE289858, has an accessible physical-DNA mapping table; its full eligibility audit remains. Keep it unused for model selection so it can provide a fresh test.

More training cycles or more AlphaGenome predictions alone are not evidence of improvement. Human genetic replication and a genuine sleep/exhaustion-excluded frailty GWAS remain separate requirements. Original genetic, NHANES and previous-model results are preserved.
'''
    (OUT/'PLAIN_LANGUAGE_RESULTS.txt').write_text(report)
    ready = json.loads((OUT/'readiness_assessment.json').read_text())
    ready['primary_molecular_test_historical'] = primary['AVI_comparison']
    ready['molecular_test'] = {'final_combined_decision': combined, 'interpretation': assessment['AVI_interpretation']}
    ready['exploratory'] = ['All unvalidated raw model projections and gene/cell annotations', 'Postselected regional genetic and cell enrichment findings', 'Posthoc zero-between-cell-difference and input-domain diagnostics']
    ready['scientific_model_success'] = assessment
    ready['next_model_plan'] = plan
    write_json(OUT/'readiness_assessment.json', ready)
    card = json.loads((OUT/'model_card.json').read_text())
    card['status'] = 'Outside tests do not establish useful magnitude prediction or reliable cell-specific ranking. All cell groups withheld.'
    card['scientific_model_success'] = assessment
    card['next_model_plan'] = plan
    write_json(OUT/'model_card.json', card)
    manuscript = ROOT/'outputs/nhanes_manuscript_audit/manuscript_with_nhanes.md'
    snapshot = OUT/'manuscript_before_V7.md'
    marker = '## Exploratory prediction of laboratory DNA effect magnitude'
    original = manuscript.read_text().split(marker)[0].rstrip()
    assert original == snapshot.read_text().rstrip(), 'Original manuscript prefix changed'
    section = '''

## Exploratory prediction of laboratory DNA effect magnitude

A separate molecular benchmark used 633,076 development observations, including 605,885 absolute log2 reporter-effect measurements and 27,191 auxiliary saturation coefficients with a separate output unit. Genomic neighborhoods and shared central sequence contexts were kept within partitions; neighborhoods overlapping the outside tests were excluded. Small neural regression heads were trained by backpropagation, with linear/tree baselines and a generic-cell comparison. The pretrained DNA encoder remained fixed. A documented development-only amendment added direct absolute-error training before outside model performance was calculated.

The primary outside test used 9,602 original fitted-effect measurements across primary neural progenitors and HEK293 cells, representing 4,820 verified physical DNA pairs from [McAfee et al. (2023)](https://doi.org/10.1016/j.xgen.2023.100404). Sequence-only mean absolute error was 0.0764 versus 0.0763 log2 units for the training-cell median in neural progenitors, and 0.1013 versus 0.0976 in HEK293 cells. Effect-magnitude Spearman correlations were −0.015 and 0.087, respectively. The additional fixed-model test used 547 baseline neuronal effects from [Retallick-Townsley et al. (2026)](https://doi.org/10.1038/s41467-026-77114-x), after duplicate reconciliation, physical-sequence verification and genomic-overlap exclusion. Error was 0.0829 versus 0.0773 for the training-cell median, with Spearman correlation 0.049. The paired genomic-window bootstrap confidence interval for error improvement over this median was entirely below zero. Thus useful cross-study magnitude prediction was not established. Original raw barcode QC and fitted-effect estimation were not independently reconstructed.

AVI, a permissive composite containing AlphaGenome-derived and other annotations, reduced macro error relative to sequence-only inputs by 2.71% in the primary study and 1.06% in the additional neuron study. The additional study failed the predeclared tolerance for reduced within-region ranking correlation; that statistic used only four eligible groups and does not establish biological harm. Both study-level gates were required, so AVI was not retained. The original primary result, which passed its individual gate, remains unchanged. These comparisons cannot attribute benefit uniquely to AlphaGenome.

Among 4,782 DNA pairs tested in both primary-study cell types, sequence-only cell-difference error was 0.1029, with correlation 0.018. A clearly labeled posthoc comparison found error 0.0909 when predicting no between-cell difference, versus 0.1029 and 0.1000 for sequence-only and AVI-added models. This comparison does not support learned variant-specific cell interactions and did not change the model or retention rule.

Nominal 90% intervals covered 97–99% of outside observations but were broad. The frozen input policy accepted none of the primary measurements because their cell groups lacked multi-study development support, and all 547 additional baseline neuronal measurements. Subsequent predeclared outside-performance protection withheld all cell groups. Coverage and abstention therefore do not establish useful precision or a universal coverage guarantee. Unsupported immune-stimulation conditions receive no authorized prediction; baseline projections compared with stimulated outcomes are diagnostic only.

The complete supplement preserves the original 1,351 gene-evidence rows and covers 498 distinct variants across five modeled cell contexts: 488 existing molecular candidates combined with 33 score markers, with 23 shared variants. There are 470 supported single-substitution inputs and 28 unsupported indels. Authorized prediction and rank fields are blank; separately labeled raw projections are unvalidated diagnostics. Nearby molecular credible-set comparisons remain unanswered. This work supplies neither person-level sleep/frailty prediction nor independent human genetic-score validation, genetic correlation, disease causation or causal target-gene evidence. The NHANES result remains separate observational evidence. [Outside-study error, coverage and answer availability](../sleep_frailty_ai_v7/figures/all_studies_error_and_answer_availability.png) are reported together; full source audits and adverse results are retained in the molecular supplement.
'''
    manuscript.write_text(original+section)
    receipt = json.loads((OUT/'manuscript_integration_receipt.json').read_text())
    receipt.update(after_sha256=sha(manuscript), original_before_section_preserved=True, original_NHANES_and_genetics_preserved=True)
    write_json(OUT/'manuscript_integration_receipt.json', receipt)
    assert before == {n: sha(OUT/n) for n in immutable}
    write_json(OUT/'scientific_report_refinement_receipt.json', {'frozen_files_unchanged': before, 'original_manuscript_prefix_preserved': True, 'new_fitting_or_threshold_change': False, 'previous_models_intact': assert_previous_intact()})
    print('REPORTS UPDATED; ADVERSE TEST RESULTS AND ORIGINAL MANUSCRIPT PRESERVED', flush=True)


if __name__ == '__main__':
    main()

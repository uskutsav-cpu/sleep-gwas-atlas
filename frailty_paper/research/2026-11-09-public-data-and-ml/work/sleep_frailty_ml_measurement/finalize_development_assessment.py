"""Record the finite negative result without consuming a new test study."""
from measurement_common import *


def main():
    assert_prior_versions()
    assessment=json.loads((DEST/'finite_measurement_assessment.json').read_text())
    numerical=json.loads((DEST/'fitted_measurement_independent_audit.json').read_text())
    counts=json.loads((DEST/'Ulirsch_count_reanalysis_audit.json').read_text())
    assert assessment['all_16_fits_complete'] and numerical['passed'] and counts['passed']
    selected=assessment['selected']
    assert selected['mean_source_line_median_error_skill'] < 0
    manifest={}
    for path in sorted((DEST/'development').glob('*_record.json')):
        r=json.loads(path.read_text()); assert r['completed_epochs']==8
        assert not r['new_Ulirsch_outcomes_used'] and not r['source_SD_used_at_inference']
        checkpoint=path.with_name(path.name.replace('_record.json','.safetensors'))
        assert sha(checkpoint)==r['checkpoint_sha256']
        for p in [path,checkpoint]: manifest[str(p.relative_to(DEST))]=sha(p)
    assert len(manifest)==32
    report={
       'completed':['16 matched development fits, two recipes, seven leave-study-out roles and one held-region role',
                    'Independent NumPy forward and objective replay; fitted finite-difference gradient audit',
                    'Fitted allele-swap and identical-DNA invariants',
                    'Original public count reconstruction and independent base-R agreement'],
       'selected_development_setting':selected,
       'decision':'Do not promote or refit this version for project predictions; it fails the primary average held-study baseline comparison.',
       'measurement_consistency_best_mean_source_skill':assessment['best_by_recipe']['measurement_consistency']['mean_source_line_median_error_skill'],
       'development_selection_is_not_independent_validation':True,
       'no_calibrated_confidence_or_patient_prediction_established':True,
       'new_count_data':{'eligible_designs':counts['eligible_designs'],'unique_variants':counts['eligible_distinct_variants'],
                         'outcomes_used_for_model_selection':False,'effect_estimates_are_new_reanalysis_not_original_author_coefficients':True,
                         'limitation':'K562 blood-cancer cell line; neither sleep/frailty outcomes nor a validation of a human genetic score.'},
       'pretrained_reader_original_benchmark_review':{
          'paper':'https://pmc.ncbi.nlm.nih.gov/articles/PMC11903340/',
          'main_text_checked':True,'Ulirsch_named_in_main_text':False,
          'reported_variant_benchmarks':'UDACHA/ADASTRA allelic imbalance and Kircher saturation mutagenesis',
          'supplementary_methods_full_review_complete':False,
          'absence_of_all_prior_benchmark_exposure_established':False},
       'next_priorities':[
          'Diagnose study-to-study sign, unit, construct, promoter and treatment differences using development studies only; retain exact known and unknown labels.',
          'Use repeat measurements to estimate which observed contrasts are reliable, without filtering test rows by significance or the desired result.',
          'Represent actual cell state and experimental conditions explicitly; neural progenitor, neuron, basal and stimulated states are not interchangeable.',
          'Before another fit, define one specific diagnostic-supported change and a finite baseline-controlled comparison; do not extend epochs or architecture searches by default.',
          'Freeze any justified candidate, preprocessing, baseline and confidence rule before opening the reserved blood study.',
          'Treat AlphaGenome outputs only as optional inputs; retain only if they improve real measured outcomes on a separate eligible study.',
          'Require whole-study and whole-region independence, useful prediction availability and relevant-cell validation before sleep/frailty mechanism claims.'
       ],
       'research_goal_complete':False,'reserved_new_study_evaluated':False,
       'genetic_and_NHANES_analyses_changed':False,
       'practical_constraint':'Internal and external disks have under 150 MiB and 100 MiB free respectively; use small metadata and existing banks. Large downloads or feature-bank expansion require storage capacity.',
       'files':manifest}
    write_json(DEST/'progress_and_readiness.json',report)
    plain='''The latest experimental model did not improve prediction across studies.

Sixteen matched development fits finished. At the selected common setting,
the best recipe made about 0.51% more absolute prediction error than the
training-derived typical effect for each cell line, averaged across seven
held-out studies. The recipe adding measurement-consistency penalties was
about 4.72% worse at its own selected common setting. These are development
comparisons, not new independent validation.

Independent numerical checks passed. That means the program calculates and
learns as specified; it does not mean its scientific predictions are useful.
This version will not replace the earlier models or generate project rankings.

An additional public blood-cell study now has 1,314 eligible DNA comparisons
covering 447 distinct variants. A separate R calculation reproduced the new
count-based estimates. These estimates remain out of model selection. They
are not the original authors' fitted effects or patient-level sleep/frailty
outcomes. Original pretrained-model supplementary benchmark exposure still
requires review before calling the entire pipeline fully unexposed.

The strongest next work is to resolve differences in experimental measurement
and actual cell state, learn from repeat experiments, and define one justified
change before more fitting. More epochs and more AI-generated target labels
are not supported by these results. AlphaGenome can be an optional input;
experimental measurements must remain the answer used to judge it.

Confidence across new experiments, real cell-specific effect differences,
independent genetic-score validation and a useful sleep-frailty prediction
system remain unfinished. Original genetic and NHANES analyses are preserved.
'''
    (DEST/'PLAIN_LANGUAGE_PROGRESS.txt').write_text(plain)
    assert_prior_versions()
    files={str(p.relative_to(DEST)):sha(p) for p in [DEST/'progress_and_readiness.json',DEST/'PLAIN_LANGUAGE_PROGRESS.txt',DEST/'finite_measurement_assessment.json',DEST/'fitted_measurement_independent_audit.json',DEST/'Ulirsch_count_reanalysis_audit.json']}
    write_json(DEST/'final_progress_verification.json',{'passed':True,'all_16_saved_weights_and_records_verified':True,
       'prior_frozen_versions_unchanged':True,'reserved_outcomes_not_evaluated':True,'files':files,
       'scientific_goal_complete':False})
    print('FINITE NEGATIVE RESULT AND NEXT PRIORITIES RECORDED',flush=True)


if __name__=='__main__': main()

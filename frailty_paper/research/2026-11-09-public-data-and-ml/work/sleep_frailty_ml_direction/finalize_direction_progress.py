"""Record the completed ablation and the unchanged scientific requirements."""
from direction_common import *


def main():
    assert_versions();a=json.loads((DEST/'finite_direction_assessment.json').read_text());s=a['selected']
    recovery=json.loads((DEST/'original_paired_signed_recovery_receipt.json').read_text())
    numerical=json.loads((DEST/'fitted_direction_independent_audit.json').read_text())
    assert a['all_seven_new_fits_complete'] and numerical['passed'] and recovery['passed']
    assert s['candidate_mean_source_skill']<0 and s['candidate_mean_source_skill']<s['control_mean_source_skill']
    epoch=s['common_epoch']-1
    curve=json.loads((DEST/'development/source_signed_none_42_curve.json').read_text())[epoch]
    old=json.loads((MEASUREMENT/'development/observation_only_none_42_curve.json').read_text())[epoch]
    manifest={}
    for p in (DEST/'development').glob('*_record.json'):
        r=json.loads(p.read_text());checkpoint=p.with_name(p.name.replace('_record.json','.safetensors'))
        assert sha(checkpoint)==r['checkpoint_sha256'] and r['completed_epochs']==8
        assert not r['new_Ulirsch_outcomes_used'] and not r['source_SD_used_at_inference']
        for item in [p,checkpoint]:manifest[str(item.relative_to(DEST))]=sha(item)
    assert len(manifest)==14
    report={
      'previous_goal_turn_classification':'Progress: completed prior finite assessment, independent fitted audit, and protected count reconstruction.',
      'current_turn_completed':['Traced original paired GLM contrast through the published source software and unambiguous construct identities',
          'Independently checked paired-design contrast algebra in base R',
          'Recovered a separate sidecar of 311006 original signed development measurements; primary targets and partitions unchanged',
          'Independent base-R exact source-coefficient check on 1799 rows across all seven original files',
          'Completed seven additional matched 8-epoch fits; wholly omitted Siraj role is exactly equivalent to the previous control',
          'Independent fitted forward, objective and gradient audit using actual newly recovered signed examples',
          'Checked a newly found cell-matched AlphaGenome-derived resource; ordinary API-derived features excluded from training under current project use terms',
          'Freed regenerable source-backed Python bytecode caches; original source files and frozen models preserved'],
      'recovered_targets':{'development':recovery['new_original_signed_development_targets'],'training':recovery['new_signed_training_rows'],
          'validation':recovery['new_signed_validation_rows'],'calibration':0,'ambiguous_excluded':66,
          'unique_new_variants_or_donors_implied':False,'historical_original_GLM_refit_established':False},
      'selected_common_epoch':s['common_epoch'],'selected_whole_study_comparison':s,
      'decision':'Reject this candidate for promotion. It does not beat the typical-effect baseline or the matched control on the primary mean held-study error endpoint.',
      'same_epoch_held_region_comparisons':{'candidate':curve['comparisons'],'control':old['comparisons']},
      'nearby_ordering_is_not_a_substitute_for_primary_effect_size_prediction':True,
      'weak_cell_specific_signal_remains':True,'new_reserved_study_outcomes_evaluated':False,
      'required_evidence':{
          'actual_effect_size_prediction':'Unmet: candidate average held-study error exceeds the baseline.',
          'cell_and_nearby_variant_differences':'Unmet: nearby ordering remains about 54%; within fixed study/cell-pair correlation is about 0.046 in development.',
          'useful_confidence_across_experiments':'Unmet: no improved outside prediction or independently useful calibration/answer availability.',
          'AlphaGenome_added_input':'Previous permissible AVI comparisons completed and rejected; the new ordinary Atlas-derived resource is not eligible for training in this project. No retained added-input improvement established.',
          'independent_genetic_score_validation':'Unmet; laboratory measurements cannot establish patient-level prediction or shared sleep-frailty genetics.'},
      'next_action_supported_by_current_evidence':'Before another full fit, test whether the existing fixed DNA representations contain learnable signed-effect signal on genuinely held regions and held studies. Compare a strongly regularized direct signed-effect readout with the existing representation; preserve actual magnitude and cell/nearby endpoints. Any representation adaptation must follow a fixed finite plan, not repeated architecture or epoch searches.',
      'reserved_source_policy':'Keep Ulirsch out of model selection. Review remaining original pretrained supplementary benchmark exposure before claiming an entirely unexposed pipeline. Freeze a justified candidate and comparison rules before any reserved-outcome evaluation.',
      'manuscript_and_observational_genetic_analyses_changed':False,
      'scientific_goal_complete':False,'files':manifest}
    write_json(DEST/'progress_and_readiness.json',report)
    (DEST/'PLAIN_LANGUAGE_PROGRESS.txt').write_text('''We recovered more real laboratory information, but the new model is still not strong enough.

The published laboratory software defines the effect as ALT-versus-REF RNA
enrichment after accounting for the corresponding DNA difference. We checked
that direction independently with base R and matched 1,799 sampled original
coefficients across all seven input files exactly. This is a direction and
source-identity audit, not an independent refit of the original GLM.

A separate file now adds 311,006 original signed development measurements,
including 232,285 training examples and 78,721 validation examples. No
calibration or new reserved-study targets were added. These are extra labels
for existing observations, not that many new DNA variants or people. Sixty-six
ambiguous source-to-construct allele mappings were excluded.

Seven additional matched training runs finished. The best common setting made
about 1.45% more absolute error than a typical-effect guess, averaged across
seven held-out studies. Its matched control at that setting was about 0.97%
worse than the guess. Nearby ordering was about 53.9%; the small improvement
there does not override failure on effect size. Cell-difference prediction
remains weak. The candidate will not generate project predictions.

Independent numerical and gradient checks passed. Those checks establish
correct calculation and actual learning, not useful scientific prediction.

We found a public cell-matched AlphaGenome-derived resource, but its ordinary
API-derived features are outside the designated permissive AVI-score training
exception. It was not used for training. The previous permissible AVI tests
and their adverse results remain unchanged.

The next diagnostic should establish whether the fixed DNA representations
contain learnable effect information before another full training exercise.
The new public blood-cell test remains out of model selection. Useful effect
prediction, real cell differences, reliable confidence and independent human
genetic-score validation are still unfinished. Original manuscript, genetic
and NHANES analyses are preserved.
''')
    # Freeze this finite result and its actual provenance, rather than modifying
    # any earlier outcome lock or claiming another weak comparison completed the goal.
    for p in sorted(DEST.glob('*.json')):
        if p.name!='final_progress_verification.json':manifest[p.name]=sha(p)
    for p in sorted(DEST.glob('*.tsv')):manifest[p.name]=sha(p)
    for p in sorted((DEST/'sources').glob('*')):manifest[str(p.relative_to(DEST))]=sha(p)
    for p in [DEST/'original_paired_signed_development.npy',DEST/'original_paired_signed_development.csv.gz',DEST/'PLAIN_LANGUAGE_PROGRESS.txt']:
        manifest[str(p.relative_to(DEST))]=sha(p)
    write_json(DEST/'final_progress_verification.json',{'passed':True,'previous_frozen_versions_unchanged':True,
       'all_seven_new_saved_weights_and_records_verified':True,'new_reserved_outcomes_not_evaluated':True,
       'candidate_promoted':False,'scientific_goal_complete':False,'files':manifest})
    assert_versions();print('SOURCE-DIRECTION ABLATION NEGATIVE RESULT VERIFIED AND PRESERVED',flush=True)


if __name__=='__main__': main()

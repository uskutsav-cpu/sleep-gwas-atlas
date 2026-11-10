"""Freeze the finite representation diagnostic and verified neuronal labels."""
from representation_common import *


def main():
    assert_versions_all();a=json.loads((DEST/'representation_diagnostic_assessment.json').read_text())
    bounds=json.loads((DEST/'convex_optimality_audit.json').read_text());assert a['all_eight_roles_complete']
    assert len(bounds['completed_roles'])==8 and a['mean_seven_held_source_magnitude_skill']<0
    rcheck=list(csv.DictReader((DEST/'independent_R_readout_numerical_receipt.tsv').open(),delimiter='\t'))[0]
    expected=json.loads((DEST/'readout_audit_expectations.json').read_text());assert rcheck['passed']=='TRUE'
    assert abs(float(rcheck['magnitude_batch_objective'])-expected['magnitude_batch_objective'])<1e-12
    folder=DEST/'neural_source_trace';nr=json.loads((folder/'neural_original_direction_recovery_receipt.json').read_text())
    nc=list(csv.DictReader((folder/'neural_direction_independent_R_receipt.tsv').open(),delimiter='\t'))[0]
    assert nc['passed']=='TRUE' and int(nc['all_recovered_rows_checked'])==nr['recovered_original_signed_targets']
    assert sha(folder/'neural_original_signed_development.npy')==nr['sidecar_sha256']
    assert sha(folder/'neural_original_signed_development.csv.gz')==nr['CSV_sha256']
    records={s:json.loads((DEST/(s+'_diagnostic.json')).read_text()) for s in get_plan()['roles']}
    for s,r in records.items(): assert sha(DEST/(s+'_readouts.npz'))==r['readouts_sha256']
    gap=max(v['fraction_of_primal'] for v in bounds['completed_roles'].values())
    write_json(DEST/'neuronal_cell_identity_review.json',{
       'source':'https://pmc.ncbi.nlm.nih.gov/articles/PMC13405209/',
       'methods_Cell_Culture_and_lentiviral_assay_reviewed':True,
       'Salomon_assay_genetic_background':'WTC11 with inducible NGN2',
       'assay_state':'Differentiated post-mitotic excitatory neurons; infected on differentiation day 7 and harvested seven days later.',
       'existing_NGN2_induced_neurons_to_WTC11_NGN2_alias_supported':True,
       'pretrained_WTC11_pluripotent_cell_reader_is_same_assay_state':False,
       'same_cell_line_name_does_not_prove_same_maturation_or_experiment':True,
       'no_training_metadata_changed':True})
    report={
       'previous_goal_turn_classification':'Progress: original paired signed-label recovery and seven real matched fits completed; adverse results preserved.',
       'current_turn_completed':['Eight fixed diagnostic roles, each with a smoothed absolute-error magnitude readout and a separate closed-form signed ridge readout',
          'Independent base-R replay of actual-data objective, gradients, finite differences and signed normal-equation solve',
          'Primal-dual optimality bounds for all eight declared convex training objectives',
          'Exact original-author Git blob verification of neuronal fitted coefficients, REF/ALT map and source CLI',
          'Recovered 29485 original neuronal signed development targets with full source construct/DNA identity checks',
          'Independent base-R direction algebra and verification of every recovered neuronal coefficient',
          'Verified neuronal genetic-background alias against original culture methods while preserving differentiation state distinction'],
       'primary_magnitude_result':a,
       'decision':'Do not promote this diagnostic readout; mean held-study magnitude error fails the training-only typical-effect baseline.',
       'optimization_limit_preserved':{'all_reported_success':a['all_magnitude_optimizations_converged'],
          'iteration_cap':80,'maximum_certified_training_objective_fractional_shortfall':gap,
          'meaning':'The finite fits hit their iteration cap. The dual bound places each within this fraction of its declared smooth penalized training objective optimum. It does not establish optimal true MAE or unseen-study accuracy.'},
       'held_region_nearby_and_cell_differences':records['none']['primary_magnitude_metrics']['comparisons'],
       'held_region_signed_diagnostics':records['none']['signed_diagnostic'],
       'interpretation':'Fixed features contain some learnable direction signal in the existing Siraj blood/liver/tumor measurements. Current neuronal direction transfer is weak. This restricted readout test does not exclude richer or nonlinear information.',
       'new_neuronal_original_signed_labels':{'total':nr['recovered_original_signed_targets'],'training':nr['training_rows'],
          'validation':nr['validation_rows'],'calibration':0,'all_independently_checked':True,'used_in_completed_diagnostic_fits':False,
          'new_independent_variants_people_or_cohorts_implied':False,'original_count_fits_reconstructed':False},
       'required_scientific_outcomes':{
          'actual_effect_size_prediction':'Unmet: average held-study error is about 0.58% worse than typical-effect guesses.',
          'cell_and_nearby_effect_differences':'Unmet: held-region nearby ordering about 54.6%; within fixed study/cell-pair correlation about 0.043.',
          'confidence_across_experiments':'Unmet: no improved external predictor or useful new confidence/answer-availability validation.',
          'AlphaGenome_added_input':'Earlier permissible AVI comparisons rejected; ordinary API-derived outputs are not new training inputs. No retained improvement established.'},
       'next_action':'Use the verified original neuronal direction labels in a finite matched neural-focused representation-adaptation comparison. Keep measured magnitude primary, actual source DNA lengths, whole-study/component exclusions and original negative results. Account for the earlier V8 two-epoch last-convolution experiment; do not claim core adaptation has never been tried. A new comparison must isolate adaptation from the new labels and protect the reserved study.',
       'no_new_test_evaluation_or_confidence_claim_or_project_ranking':True,
       'original_manuscript_NHANES_and_genetic_results_changed':False,'scientific_goal_complete':False}
    write_json(DEST/'progress_and_readiness.json',report)
    (DEST/'PLAIN_LANGUAGE_PROGRESS.txt').write_text('''We completed a check of whether the existing DNA features carry useful effect information.

Eight fixed development roles finished. Each role tested a simple effect-size
readout and a separate direction readout. Averaged across seven held-out
studies, size error was about 0.58% worse than a training-derived typical-effect
guess. The model is not promoted. Nearby ordering was about 54.6% on held
regions; fine cell-difference prediction remains weak.

The simple fits reached their fixed 80-iteration limit. An independent
mathematical bound shows each is within 0.014% of the best value of its chosen
smooth, penalized training objective. That does not prove minimum true error
or successful prediction in another study. A separate base-R check reproduced
the calculations and checked the actual-data gradients and ridge solution.

There is some learnable direction signal in existing blood/liver/tumor-cell
measurements. The neuronal direction signal in these features is much weaker.
This does not prove all possible features or nonlinear models would fail.

We also recovered 29,485 original neuronal effect-direction measurements for
existing observations: 21,577 training and 7,908 validation rows. Every recovered
coefficient passed an independent R comparison with the original fitted-results
table. Original REF/ALT construct identifiers and entire physical DNA pairs
matched our inputs. No calibration labels were added. These new labels were
not used in the just-completed diagnostic fits.

The source culture methods confirm WTC11-derived NGN2 neurons. The existing
cell-line alias is supported, while a WTC11 stem-cell model remains a different
cell state. We recovered original fitted coefficients, not independent raw
count fits or the exact historical software environment.

The next justified comparison is a focused adaptation using these verified
neuronal labels, with the same labels in both adaptation and frozen controls.
Effect size stays the primary outcome. Earlier core-tuning attempts and all
negative results must remain visible. The reserved new study stays out of
model selection. No useful disease-risk model, causal mechanism, independent
genetic-score validation or completed scientific goal is established.
''')
    manifest={}
    for p in sorted(DEST.rglob('*')):
        if p.is_file() and p.name!='final_progress_verification.json':manifest[str(p.relative_to(DEST))]=sha(p)
    code={p.name:sha(p) for p in Path(__file__).parent.glob('*.py')}
    code.update({p.name:sha(p) for p in Path(__file__).parent.glob('*.R')})
    write_json(DEST/'final_progress_verification.json',{'passed':True,'files':manifest,'source_code_sha256':code,
       'previous_frozen_versions_unchanged':True,'all_eight_readout_artifacts_verified':True,
       'new_neuronal_labels_not_used_in_completed_fits':True,'reserved_new_study_not_evaluated':True,
       'candidate_promoted':False,'scientific_goal_complete':False})
    assert_versions_all();print('FINITE REPRESENTATION RESULT AND NEW NEURONAL LABELS VERIFIED',flush=True)


if __name__=='__main__':main()

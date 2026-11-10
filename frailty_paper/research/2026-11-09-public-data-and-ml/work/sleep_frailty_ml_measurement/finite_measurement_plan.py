"""One bounded comparison of explicit observed-effect and measurement consistency."""
from measurement_common import *

def get_plan():
    plan={'data_sha256':sha(OLD/'coverage_training_arrays.npz'),'pairs_sha256':sha(OLD/'coverage_measured_pairs.npz'),'channel_receipt_sha256':sha(CHANNEL/'channel_feature_receipt.json'),
          'primary_endpoint':'Original measured absolute log2 effect magnitude. Same context-balanced MAE as prior versions, unmodified by SE; detection is not a training target. No denoised component replaces this endpoint.',
          'architecture':'Shared cell/line-conditioned odd signed-effect branch over fixed per-reader global/local mean and difference projections, symmetric mean-activity branch and positive assay-noise branch. Observed magnitude is the median of abs(N(predicted signed effect,predicted assay SD)).',
          'structural_invariants':'Allele swapping negates the signed effect and preserves magnitude; identical allele inputs have zero signed component. The observed-magnitude component can remain positive due to assay variability.',
          'recipes':['observation_only','measurement_consistency'],'held_studies':['none','Salomon','Kosicki2025','Siraj','IGVF240K','Koesterich2023','McAfee2023','Retallick2026'],'epochs':8,'seed':42,
          'comparison':'Identical architecture, initial parameters, rows, normalization, batches, primary/RNA/paired objectives and update count. Only measurement_consistency adds the two declared measurement auxiliaries.',
          'extra_auxiliaries':'0.05 robust folded-Gaussian consistency with original measured SE, and 0.05 smooth-L1 agreement of log predicted assay SD with log source SE. Apply context/precision weights as other auxiliaries; original main magnitude MAE remains unchanged. Missing SE stays absent.',
          'auxiliary_assumptions':'Reported SE approximates coefficient sampling variation; Gaussian folded observation is a working approximation, with square-root robust deviance for large residuals. Source SE floor 0.001 is numerical regularization. Neither component is validated biological truth or a calibrated confidence interval.',
          'normalization':'Mean-channel location and SD estimated from unique training physical banks only; signed-difference channels use training RMS without centering to preserve allele antisymmetry. Physical length/position normalization uses training banks only.',
          'source_transforms':'Same earlier training-only source/cell activity centers and SE precision scaling. No current measured RNA, SE, p value, study ID, experiment ID or new-study label is supplied at inference.',
          'selection':'One common recipe/epoch maximizing mean seven whole-source-omission measured-magnitude MAE improvement over training-only cell/culture median. Report local and within-fixed-cell-pair diagnostics separately; no favorable per-source epochs.',
          'stopping_rule':'Complete exactly this comparison; no epoch, source or architecture expansion inside it. Do not promote from a better auxiliary, latent, detection or within-source score alone.',
          'new_test_policy':'Ulirsch outcome coefficients are sealed and may be joined only after a justified candidate, confidence rules and predictions are fixed. Previous blood/OL33/Matoba tests are already revealed for any changed model.',
          'prior_AVI_rejection':'Preserved; no AVI or restricted AlphaGenome teacher values are used in this development comparison.',
          'pretrained_reader_weights_updated':False,'scientific_objective_complete':False}
    path=DEST/'finite_measurement_plan.json'
    if path.exists():assert json.loads(path.read_text())==plan
    else:write_json(path,plan)
    return plan

if __name__=='__main__':
    assert_prior_versions();get_plan();print('FINITE MEASUREMENT-AWARE PLAN FIXED',flush=True)

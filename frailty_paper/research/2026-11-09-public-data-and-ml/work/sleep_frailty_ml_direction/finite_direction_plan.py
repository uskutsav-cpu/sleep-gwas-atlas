from direction_common import *


def get_plan():
    plan={'recipes':['source_signed'],'held_studies':['none','Salomon','Kosicki2025','Siraj','IGVF240K','Koesterich2023','McAfee2023','Retallick2026'],
       'epochs':8,'seed':42,'data_sha256':sha(OLD/'coverage_training_arrays.npz'),
       'signed_sidecar_sha256':sha(DEST/'original_paired_signed_development.npy'),
       'change':'Only populate documented original paired GLM signed targets for existing eligible Siraj training/development observations.',
       'primary_endpoint':'Actual observed magnitude MAE versus training-only typical cell-line effect, macro across held studies; no SE reweighting of primary error.',
       'architecture_and_optimization':'Same FoldedEffectHead, inputs, 8 epochs, initial parameters, batches, RNA/local/cross objectives, optimizer and seed as observation_only control. No measurement-consistency penalty.',
       'control':'Existing observation_only development fits in sleep_frailty_ml_measurement; immutable and numerically audited.',
       'Siraj_omission':'Reuse its existing observation_only curve because every new signed label belongs to the wholly omitted study; no additional fitting can change this training role.',
       'selection':'One common epoch using mean seven whole-study-omission magnitude skill. Report matched control at each common epoch and both complete curves.',
       'stopping_rule':'Exactly this one supervision comparison. Do not extend epochs, search loss weights, or open Ulirsch to rescue a negative result.',
       'reserved_outcomes_used':False,'calibration_labels_changed':False,'AlphaGenome_outputs_added':False,
       'original_frozen_models_and_primary_targets_changed':False,
       'interpretation':'Direction recovery and a successful fit do not establish external predictive improvement, causal genes, confidence calibration or human risk prediction.'}
    path=DEST/'finite_direction_plan.json'
    if path.exists(): assert json.loads(path.read_text())==plan
    else: write_json(path,plan)
    return plan


if __name__=='__main__': assert_versions();get_plan();print('FINITE ORIGINAL-DIRECTION COMPARISON PLAN FIXED')

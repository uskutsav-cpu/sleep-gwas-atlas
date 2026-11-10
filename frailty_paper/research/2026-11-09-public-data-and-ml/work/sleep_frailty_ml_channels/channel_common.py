"""Separate development version; all prior frozen results remain immutable."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_next'))
OLD=ROOT/'outputs/sleep_frailty_ml_next'
DEST=ROOT/'outputs/sleep_frailty_ml_channels'
DEST.mkdir(exist_ok=True)

def assert_frozen_previous():
    assert_previous_sources()
    lock=OLD/'validation/before_reserved_outcomes_lock.json'
    for path,digest in json.loads(lock.read_text())['files'].items():assert sha(path)==digest,path
    result=OLD/'validation/reserved_outcome_evaluation.json'
    assert sha(result)==json.loads((OLD/'validation/reserved_results_independent_audit.json').read_text())['original_result_sha256']

def get_plan():
    plan={'data':'Unchanged 654,612 observations in prior coverage_training_arrays.npz; no blood or OL33 outcomes added.',
          'source_data_sha256':sha(OLD/'coverage_training_arrays.npz'),'source_pairs_sha256':sha(OLD/'coverage_measured_pairs.npz'),
          'features':'Same 66 scalar/physical inputs plus 16 fixed orthogonal projections of each reader global/local mean and signed allele difference. 192 additional coordinates from three frozen 128-channel readers.',
          'projection_dimension':16,'projection_seed':20261108,'difference_storage_multiplier':1000.,'storage_dtype':'float16; cast to float32 and normalize using training-only physical banks.',
          'recipes':['scalar_mask','channel_preserving'],'held_studies':['none','Salomon','Kosicki2025','Siraj','IGVF240K','Koesterich2023','McAfee2023','Retallick2026'],
          'epochs':8,'seed':42,'architecture':'Same CellMatchedEffectHead with 258 inputs in both arms. Scalar control masks 192 added columns; identical head initialization, data, loss, batches, dropout schedule and updates.',
          'objective':'Primary actual magnitude MAE unchanged and unmodified by SE; same auxiliary activity and signed targets and paired nearby/cell-difference losses. No detected/not-detected target.',
          'selection':'One common recipe/epoch by mean seven-whole-source-omission MAE skill over training-only cell/culture medians; held-region pairs reported separately; no per-source selected epoch presented as independent validation.',
          'stopping_rule':'Complete this one finite information-preservation comparison. Do not expand epochs, projection dimensions or sources within it. No final ensemble/confidence refit or fresh-study test without evidence justifying the selected candidate and a newly specified validation plan.',
          'prior_tests':'Blood, OL33, Matoba, McAfee and Retallick previously observed results remain at their original frozen versions; none is untouched for a changed model.',
          'new_untouched_test':'Not yet established. Independently audit original constructs, coefficients, source history and pretrained exposure before claiming new-study validation.',
          'AVI':'Previously rejected. Not used in this development comparison; no claim to have fine-tuned AlphaGenome.',
          'pretrained_reader_weights_updated':False,'hypothesis':'Scalar summaries may discard useful channel identity. A successful implementation is not evidence that this hypothesis or scientific objective is true.'}
    path=DEST/'finite_channel_plan.json'
    if path.exists():assert json.loads(path.read_text())==plan
    else:write_json(path,plan)
    return plan

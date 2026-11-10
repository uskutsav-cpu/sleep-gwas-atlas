import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'sleep_frailty_ml_bounded'))
from bounded_common import *
BOUNDED=DEST
DEST=ROOT/'outputs/sleep_frailty_ml_resumable'
original_assert_versions=assert_versions_stable
sys.path.insert(0,str(Path(__file__).resolve().parent))

def assert_versions_stable():
    original_assert_versions()
    manifest=json.loads((BOUNDED/'interrupted_run_manifest.json').read_text())
    for path,h in manifest['files_sha256'].items():assert sha(ROOT/path)==h,path
    receipt=json.loads((DEST/'source_derivation_receipt.json').read_text())
    for path,h in receipt['previous_code_sha256'].items():assert sha(ROOT/path)==h,path

def get_stable_plan():
    plan=json.loads((BOUNDED/'finite_stable_native_plan.json').read_text())
    plan['recovery_execution']={
        'prior_plan_sha256':sha(BOUNDED/'finite_stable_native_plan.json'),
        'interruption_manifest_sha256':sha(BOUNDED/'interrupted_run_manifest.json'),
        'checkpoint_every_updates':100,
        'checkpoint_contents':'Trainable parameters, Adam moments/steps, CPU and MPS RNG, NumPy sampler state, saved epoch losses and completed curves; atomic replacement; weights_only=True loading.',
        'resumption':'Reconstruct and verify the original batch/pair schedule from seed without optimizer steps for past batches. Restore state before the next real update. Preserve artifacts beyond the last committed checkpoint before replay.',
        'initial_recovery':'Original interrupted run lacks optimizer state. This version starts again from the same pinned initial weights, preserving the partial original; it does not initialize from its partial weights or select between repeats.',
        'scientific_labels_model_objective_epochs_seed_and_gates_unchanged':True,
        'new_neuronal_activity_inventory_used_for_fit':False}
    path=DEST/'finite_stable_native_plan.json'
    if path.exists():assert json.loads(path.read_text())==plan
    else:write_json(path,plan)
    return plan

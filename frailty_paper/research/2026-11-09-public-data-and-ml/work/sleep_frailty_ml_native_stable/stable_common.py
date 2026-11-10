import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'sleep_frailty_ml_native'))
from native_common import *
ORIGINAL_NATIVE=DEST
DEST=ROOT/'outputs/sleep_frailty_ml_native_stable'


def assert_versions_stable():
    assert_versions_native()
    for path,h in json.loads((ORIGINAL_NATIVE/'resource_failure_manifest.json').read_text())['preserved_files_sha256'].items():assert sha(ROOT/path)==h,path


def get_stable_plan():
    plan=json.loads((ORIGINAL_NATIVE/'finite_native_adaptation_plan.json').read_text())
    plan['resource_revision']={
       'original_plan_sha256':sha(ORIGINAL_NATIVE/'finite_native_adaptation_plan.json'),
       'failure_manifest_sha256':sha(ORIGINAL_NATIVE/'resource_failure_manifest.json'),
       'reason':'Confirmed memory growth to 5.8G and disk exhaustion while full swap expanded; varying-shape MPS graph-cache growth is a supported hypothesis, not a proven sole cause.',
       'batching':'Shuffle actual rows within each true-length bucket, split into main batches, then shuffle the batch list. Each main training row appears exactly once per epoch. Pair samples use the same true length. No padding or invented bases.',
       'no_within_batch_construct_deduplication':'Native reader evaluates repeated actual constructs separately to keep shapes fixed; frozen BN/eval makes outputs independent of duplicate batch membership.',
       'loss':'Algebraically identical masked means using fixed-size operations; missing labels retain zero loss weight and are never assigned an observed zero effect.',
       'validation':'Predict in true-length buckets and restore the exact original observation order.',
       'device':'cpu',
       'memory_policy':'CPU execution with two threads after the bounded-shape MPS preflight exhausted its 2.13 GiB allocation limit. Preserve that failed preflight and original resource plan. Run the full 1000-step bounded-shape CPU audit before scientific fits; enforce physical footprint and disk reserve limits.',
       'original_partial_run_retained':True,'hyperparameters_epochs_targets_and_model_architecture_unchanged':True,
       'comparison_restarts_from_same_original_author_weights':True}
    path=DEST/'finite_stable_native_plan.json'
    if path.exists():assert json.loads(path.read_text())==plan
    else:write_json(path,plan)
    return plan

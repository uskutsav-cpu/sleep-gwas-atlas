import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'sleep_frailty_ml_native_stable'))
from stable_common import *
STABLE=DEST
DEST=ROOT/'outputs/sleep_frailty_ml_bounded'
original_assert_versions=assert_versions_stable

def assert_versions_stable():
    original_assert_versions()
    receipt=json.loads((DEST/'source_derivation_receipt.json').read_text())
    for path,h in receipt['preserved_previous_sources_and_attempts'].items():assert sha(ROOT/path)==h,path

def get_stable_plan():
    plan=json.loads((STABLE/'finite_stable_native_plan.json').read_text())
    plan['bounded_native_reader_execution']={
       'prior_CPU_plan_sha256':sha(STABLE/'finite_stable_native_plan.json'),
       'CPU_stop_receipt_sha256':sha(STABLE/'incomplete_CPU_preflight_stop.json'),
       'native_microbatch':32,
       'unused_rows':'Fill a short native-reader microbatch only with exact copies of its last real source construct. Remove those extra outputs before feature construction, head evaluation and losses. No DNA bases, flanks, measured targets or supervised observations are invented.',
       'equivalence':'Reader BN stays in eval mode and has no dropout; no operation mixes observations. Head is called once on the original complete actual observation batch. Discarded copies contribute no gradient.',
       'resource_policy':'MPS high/low watermark 0.4/0.3; real native reader input always 64 allele sequences by 4 bases by actual DNA length. Full 1000-update sustained preflight, all real length/tail/pair shapes, final complete validation pass.',
       'scientific_hyperparameters_loss_data_and_epochs_unchanged':True}
    plan['resource_revision']['device']='mps'
    path=DEST/'finite_stable_native_plan.json'
    if path.exists():assert json.loads(path.read_text())==plan
    else:write_json(path,plan)
    return plan

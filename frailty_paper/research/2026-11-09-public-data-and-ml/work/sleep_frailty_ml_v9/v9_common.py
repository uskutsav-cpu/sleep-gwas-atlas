"""New data-coverage experiment; all prior versions remain immutable."""
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parents[1]/'sleep_frailty_ml_v8'))
from v8_common import ROOT, BASE, V7, V6, sha, write_json, write_csv, csv_rows, rc, assert_previous_intact, environment
import os, csv, gzip, json, hashlib, collections
import numpy as np

WORK = Path(__file__).resolve().parent
V8 = BASE/'sleep_frailty_ml_v8'
DATA = BASE/'sleep_frailty_ml_v9'
OUT = DATA/'outputs'
MODELS = DATA/'models'
for p in [DATA, OUT, MODELS]: p.mkdir(parents=True, exist_ok=True)
alias = ROOT/'outputs/sleep_frailty_ai_v9'
if not alias.exists(): alias.symlink_to(OUT, target_is_directory=True)
assert alias.resolve() == OUT.resolve()


def assert_previous_sources():
    assert_previous_intact()
    r = json.loads((V8/'outputs/measured_activity_uncertainty_enrichment_audit.json').read_text())
    assert sha(V8/'measured_training_arrays_with_activity_and_SE.npz') == r['new_arrays_sha256']
    for name, receipt in [('pretrained_pooled128.npy', 'pretrained_feature_receipt.json'), ('pretrained_variant_local128.npy', 'variant_local_feature_receipt.json')]:
        assert sha(V8/name) == json.loads((V8/'outputs'/receipt).read_text())['sha256']

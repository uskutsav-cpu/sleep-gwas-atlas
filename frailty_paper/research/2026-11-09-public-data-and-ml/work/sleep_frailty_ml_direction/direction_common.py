import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'sleep_frailty_ml_measurement'))
from measurement_common import *
MEASUREMENT=DEST
MEASURE_WORK=ROOT/'work/sleep_frailty_ml_measurement'
DEST=ROOT/'outputs/sleep_frailty_ml_direction'


def assert_versions():
    assert_prior_versions()
    r=json.loads((MEASUREMENT/'final_progress_verification.json').read_text())
    for name,h in r['files'].items(): assert sha(MEASUREMENT/name)==h,name


def load_signed_sidecar(d):
    receipt=json.loads((DEST/'original_paired_signed_recovery_receipt.json').read_text())
    assert receipt['passed'] and sha(DEST/'original_paired_signed_development.npy')==receipt['target_sha256']
    assert sha(OLD/'coverage_training_arrays.npz')==receipt['original_training_arrays_sha256']
    audit=list(csv.DictReader((DEST/'signed_source_independent_R_receipt.tsv').open(),delimiter='\t'))
    assert len(audit)==len(receipt['source_files']) and all(r['passed']=='TRUE' for r in audit)
    signed=np.load(DEST/'original_paired_signed_development.npy');good=np.isfinite(signed)
    assert len(signed)==len(d['effect']) and not np.isfinite(signed[d['split']==2]).any()
    assert not (good&np.isfinite(d['signed_effect'])).any()
    assert np.allclose(abs(signed[good]),d['effect'][good],rtol=1e-5,atol=1e-6)
    d['signed_effect']=d['signed_effect'].copy();d['signed_effect'][good]=signed[good]
    return good

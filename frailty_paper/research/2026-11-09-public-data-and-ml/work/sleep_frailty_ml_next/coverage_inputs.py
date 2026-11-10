"""Read-only assembly of the frozen old bank and small appended feature bank."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
DEST=ROOT/'outputs/sleep_frailty_ml_next'

def load_inputs():
    receipt=json.loads((DEST/'coverage_expansion_receipt.json').read_text());features=json.loads((DEST/'coverage_feature_receipt.json').read_text())
    assert receipt['passed'] and features['passed']
    for path,h in [('coverage_training_arrays.npz',receipt['training_arrays_sha256']),('coverage_measured_pairs.npz',receipt['pairs_sha256']),
                   ('coverage_appended_reader_features.npy',features['features_sha256']),('coverage_appended_physical_metadata.npy',features['physical_metadata_sha256'])]:assert sha(DEST/path)==h
    old=json.loads((DEST/'cell_matched_feature_receipt.json').read_text());assert sha(DEST/'cell_matched_native_development_features.npy')==old['features_sha256'];assert sha(DEST/'native_development_physical_metadata.npy')==old['physical_metadata_sha256']
    with np.load(DEST/'coverage_training_arrays.npz') as z:d={k:z[k] for k in z.files}
    with np.load(DEST/'coverage_measured_pairs.npz') as z:pairs={k:z[k] for k in z.files}
    bank=np.concatenate([np.load(DEST/'cell_matched_native_development_features.npy',mmap_mode='r'),np.load(DEST/'coverage_appended_reader_features.npy')])
    physical=np.concatenate([np.load(DEST/'native_development_physical_metadata.npy'),np.load(DEST/'coverage_appended_physical_metadata.npy')]);used=np.unique(d['bank_index'])
    assert np.isfinite(bank[used]).all();x=np.full((len(bank),66),np.nan,dtype='float32');x[used,:48]=bank[used].reshape(len(used),48)
    pp=physical[used];assert np.all(pp[:,0]>0);x[used,48]=pp[:,0]/230.;x[used,49]=pp[:,1]/np.maximum(pp[:,0]-1,1);x[used,50:]=np.eye(16,dtype='float32')[pp[:,2]*4+pp[:,3]]
    return d,pairs,bank,x

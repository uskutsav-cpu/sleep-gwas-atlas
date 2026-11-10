"""Reload locally fitted v2 heads and comparators. Never fit on inference data."""
import pickle,numpy as np
from common import MODELS,OUT

def infer_record(record,features,records,pairs=None):
    family=record['family'];name=record['name'];feature=record['features']
    if family=='neural_seed_ensemble':return np.mean([infer_record(records[n],features,records,pairs) for n in record['config']['members']],axis=0)
    if family=='direct_avi':
        # Training median is stored once in the selected preprocessing receipt.
        import json
        med=json.loads((OUT/'preprocessing_receipt.json').read_text())['avi_training_median_raw']
        x=features['avi'][:,0];x=np.where(np.isfinite(x),x,med)
        return 1/(1+np.exp(-np.clip(x,-30,30)))
    if family=='v1_logistic':
        if pairs is None:raise ValueError('Original paired sequences required for v1 comparator')
        from common import OLD,ROOT
        import sys
        sys.path.insert(0,str(OLD));from predict import predict
        return predict(ROOT/'outputs/sleep_frailty_ai_v1/models/kmer_baseline.npz',pairs,'numpy')
    if family=='distilled_hyenadna':
        if pairs is None:raise ValueError('Paired sequences required for the native-weight DNA student')
        from common import OLD
        import sys
        sys.path.insert(0,str(OLD));from predict import predict
        return predict(OUT/record['artifact'],pairs,'mlx')
    path=OUT/record['artifact'];x=features[feature]
    if family in ['logistic','hist_gradient_boosting']:
        # Only our own local trusted checkpoint, never untrusted pickle files.
        with path.open('rb') as f:model=pickle.load(f)
        return model.predict_proba(x)[:,1]
    if family=='neural_head':
        z=np.load(path);x=np.where(np.isfinite(x),x,z['imputation']);x=np.clip((x-z['mean'])/z['sd'],-20,20).astype('float32')
        logits=np.tanh(x@z['w1']+z['b1'])@z['w2']+z['b2']
        return 1/(1+np.exp(-np.clip(logits.astype('float64'),-30,30)))
    raise ValueError(family)

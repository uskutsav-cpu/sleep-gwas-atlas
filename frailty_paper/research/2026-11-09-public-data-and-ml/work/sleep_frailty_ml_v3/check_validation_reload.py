"""Check completed candidate checkpoints against saved validation predictions.

Only chromosomes4–5 are read for inference. No blind test is unlocked here.
"""
import json,numpy as np
from common import *
from validation_registry import completed_candidates
from inference_records import StoredFeatureInference

records,paths=completed_candidates();assert records
engine=StoredFeatureInference();ix=engine.features.validation;assert set(engine.data['chromosome'][ix])=={4,5};checks={}
for name,r in records.items():
    with np.load(paths[name]) as z:
        assert np.array_equal(z['row_index'],ix);expected=z[name]
    actual=engine.predict(r,ix,records);error=float(np.max(np.abs(actual['overall']-expected)));assert error<3e-5,(name,error)
    measured=study_metrics(engine.data['label'][ix],actual['overall'],engine.data['study'][ix]);delta={s:abs(measured['studies'][s]['average_precision']-r['validation']['studies'][s]['average_precision']) for s in STUDIES}
    assert max(delta.values())<1e-4,(name,delta)
    extra={}
    if 'conditional_scores' in actual:
        prefix=paths[name].name.removesuffix('_validation_predictions.npz');cp=OUT/(prefix+'_conditional_validation_predictions.npz')
        with np.load(cp) as z:assert np.array_equal(z['row_index'],ix);conditional_expected=z[name]
        conditional_error=float(np.max(np.abs(actual['conditional_scores'][:,1]-conditional_expected)));assert conditional_error<3e-5,(name,conditional_error);extra['conditional_prediction_max_error']=conditional_error
    checks[name]={'max_absolute_prediction_error':error,'AP_absolute_difference':delta,**extra}
    if len(checks)%10==0:print('Validation checkpoint replay',len(checks),'/',len(records),flush=True)
(OUT/'validation_reload_check.json').write_text(json.dumps({'models_checked':len(checks),'scope':'Independent NumPy reload on validation chromosomes4-5; not biological validation or blind test performance','checks':checks,'blind_test_evaluated':False},indent=2)+'\n');print('VALIDATION RELOAD CHECK COMPLETE',len(checks),flush=True)

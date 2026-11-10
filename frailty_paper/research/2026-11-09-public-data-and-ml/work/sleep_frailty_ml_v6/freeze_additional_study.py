"""Freeze the extra complete study's predictions before measured-effect comparison."""
import datetime
from common import *
from encode_inputs import features
from predict_checkpoints import predict_registry
from calibrate_uncertainty import apply_policy
def main():
    assert not (OUT/'additional_test_evaluation_started.json').exists()
    original=json.loads((OUT/'final_model_lock.json').read_text())
    for path,h in original['files'].items():assert sha(path)==h,path
    for path,h in original['model_checkpoints'].items():assert sha(path)==h,path
    a=json.loads((OUT/'additional_study_blind_design_audit.json').read_text());p=DATA/'additional_study_blind_design.csv';assert sha(p)==a['design_sha256'];rows=list(csv.DictReader(p.open()));fp=DATA/'additional_study_fixed_features.npy';receipt=OUT/'additional_study_fixed_features_receipt.json'
    if fp.exists() and receipt.exists():assert sha(fp)==json.loads(receipt.read_text())['sha256']
    else:
        x=features([r['reference'] for r in rows],[r['alternate'] for r in rows],[int(r['changed_position']) for r in rows]);np.save(fp,x);write_json(receipt,{'sha256':sha(fp),'shape':list(x.shape),'design_sha256':sha(p),'measured_effects_used':False})
    x=np.load(fp,mmap_mode='r');cells=np.full(len(rows),'liver');registry=json.loads((OUT/'final_fit_registry.json').read_text());policy=json.loads((OUT/'uncertainty_policy.json').read_text());prediction,members,distance=predict_registry(registry,x,np.arange(len(rows)),cells);confidence=apply_policy(prediction['primary'],cells,distance,policy)
    pred=OUT/'additional_study_blinded_predictions.npz';np.savez_compressed(pred,name=np.array([r['name'] for r in rows]),distance=distance,**{'model::'+k:v for k,v in prediction.items()},**{'member::'+k:v for k,v in members.items()},**confidence)
    critical=[p,fp,receipt,pred,OUT/'additional_study_blind_design_audit.json',OUT/'additional_reference_audit.json',WORK/'ADDITIONAL_VALIDATION_PROTOCOL.txt',WORK/'prepare_additional_study.py',WORK/'additional_reference.py',WORK/'freeze_additional_study.py',WORK/'evaluate_additional_study.py',SOURCES/'recovery_IGVFFI7902MFBK.gz',SOURCES/'recovery_IGVFFI4642CBGA.gz',SOURCES/'recovery_IGVFFI2260ZLFH.tsv.gz',SOURCES/'recovery_IGVFFI1642PFRO.gz']
    write_json(OUT/'additional_study_validation_lock.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'original_model_lock_sha256':sha(OUT/'final_model_lock.json'),'files':{str(p):sha(p) for p in critical},'prediction_rows':len(rows),'models_thresholds_unchanged':True,'additional_study_model_outcome_comparison_already_performed':False,'qualification':'Availability-driven extra validation following broken canonical original test output. Entire original model/policy lock retained; no fitting or selection after test opening. Predictions fixed before comparison with extra study outcomes.'});print('ADDITIONAL STUDY PREDICTIONS FROZEN',len(rows),flush=True)
if __name__=='__main__':main()

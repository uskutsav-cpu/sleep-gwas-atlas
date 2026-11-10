"""Freeze models, uncertainty rules, code and outcome-blind test predictions."""
import subprocess,datetime
from common import *
from predict_checkpoints import predict_registry
from calibrate_uncertainty import apply_policy
def main():
    assert not (OUT/'final_test_evaluation_started.json').exists()
    for name in ['full_feature_gradient_check.json','final_checkpoint_audit.json','scorer_smoke_audit.json']:
        assert json.loads((OUT/name).read_text())['passed'],name
    processes=subprocess.run(['ps','ax','-o','comm=,args='],capture_output=True,text=True,check=True).stdout.splitlines()
    for line in processes:
        comm=line.strip().split(None,1)[0] if line.strip() else ''
        if comm.lower().endswith(('python','python3','python3.13')):
            assert not any('sleep_frailty_ml_v5/'+s in line for s in ['train_development.py','select_and_refit.py']),line
    registry=json.loads((OUT/'final_fit_registry.json').read_text());policy=json.loads((OUT/'uncertainty_policy.json').read_text());assert sha(OUT/'final_fit_registry.json')==policy['primary_registry_sha256'];assert {'local_mlp','contrast_standard','contrast_hard','contrast_neural_only'}<={r['name'] for r in registry['records']}
    models={}
    for r in registry['records']:assert sha(r['checkpoint'])==r['checkpoint_sha256'];models[r['checkpoint']]=r['checkpoint_sha256']
    receipt=json.loads((OUT/'kosicki_fixed_features_receipt.json').read_text());assert sha(receipt['file'])==receipt['sha256'];rows=list(csv.DictReader((DATA/'kosicki_blind_design.csv').open()));x=np.load(receipt['file'],mmap_mode='r');assert len(x)==len(rows);cells=np.full(len(rows),'neuronal_induced')
    predictions,members,distance=predict_registry(registry,x,np.arange(len(rows)),cells);uncertainty=apply_policy(predictions['primary'],cells,distance,policy);path=OUT/'kosicki_blinded_predictions.npz';np.savez_compressed(path,name=np.array([r['name'] for r in rows]),distance=distance,**{'model::'+k:v for k,v in predictions.items()},**{'member::'+k:v for k,v in members.items()},**uncertainty)
    critical=[OUT/n for n in ['final_fit_registry.json','uncertainty_policy.json','fixed_features_receipt.json','kosicki_fixed_features_receipt.json','development_selection.json','model_search_plan.json','partition_audit.json','full_feature_gradient_check.json','final_checkpoint_audit.json','scorer_smoke_audit.json','selected_motifs.json','uncertainty_protocol_details.json']]+[WORK/n for n in ['FINAL_EVALUATION_PROTOCOL.txt','EXPERIMENT_PLAN.txt','network.py','predict_checkpoints.py','score.py','encode_inputs.py','calibrate_uncertainty.py','prepare_features.py','data.py','metrics.py','train_development.py','select_and_refit.py']]+[DATA/'kosicki_blind_design.csv',SOURCES/'kosicki_mpravista.zip',SOURCES/'kosicki_sealed_outcome_manifest.json',OLD/'hyenadna/model.safetensors',OLD/'model.py',OLD/'model_mlx.py',ROOT/'work/sleep_frailty_ml_v3/motif_model.py',ROOT/'work/sleep_frailty_ml_v3/motifs/H13CORE_pwms.txt',path]
    frozen={str(p):sha(p) for p in critical};write_json(OUT/'final_model_lock.json',{'frozen_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':frozen,'model_checkpoints':models,'prediction_rows':len(rows),'final_outcomes_deserialized':False,'purpose':'Unblinding lock: all model choices, weights, confidence rules and test predictions fixed before measured outcomes are opened.','V3_registry_sha256':sha(ROOT/'outputs/sleep_frailty_ai_v3/final_model_registry.json')});print('FINAL MODEL AND OUTCOME-BLIND PREDICTIONS FROZEN',len(rows),flush=True)
if __name__=='__main__':main()

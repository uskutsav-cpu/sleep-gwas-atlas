"""Exercise the actual CSV scorer, including weak/unsupported/invalid inputs."""
import gzip,subprocess,sys
from common import *
from data import load
from predict_checkpoints import predict_registry
from calibrate_uncertainty import apply_policy
def main():
    assert not (OUT/'final_test_evaluation_started.json').exists();d=load();bank=d['bank'];chosen=[]
    for name in ['neuronal_induced','liver','colon','blood']:
        ids=np.flatnonzero((d['split']==2)&np.array([d['cell_names'][i]==name for i in d['cell']]))[:2]
        chosen.extend(ids.tolist())
    assert len(chosen)==8;rows=[];bi=d['bank_index'][chosen];cells=np.array([d['cell_names'][i] for i in d['cell'][chosen]])
    for i,index in enumerate(chosen):
        b=bank[int(d['bank_index'][index])];rows.append({'id':'calibration_input_'+str(index),'reference':b['reference'],'alternate':b['alternate'],'cell_group':str(cells[i])})
    rows.append({**rows[0],'id':'unsupported_cell_probe','cell_group':'unseen_cell_group'});rows.append({'id':'invalid_sequence_probe','reference':'N'*151,'alternate':'N'*150+'A','cell_group':'neuronal_induced'})
    src=DATA/'scorer_smoke_inputs.csv';dst=OUT/'scorer_smoke_outputs.csv'
    with src.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    subprocess.run([sys.executable,str(WORK/'score.py'),'--input',str(src),'--output',str(dst)],check=True)
    got=list(csv.DictReader(dst.open()));assert len(got)==10;assert got[-1]['status']=='unranked_invalid_input' and got[-1]['ranked_score']=='';assert got[-2]['status']=='unranked_unsupported_cell' and got[-2]['ranked_score']==''
    receipt=json.loads((OUT/'fixed_features_receipt.json').read_text());x=np.load(receipt['feature_path'],mmap_mode='r');registry=json.loads((OUT/'final_fit_registry.json').read_text());policy=json.loads((OUT/'uncertainty_policy.json').read_text());pred,members,distance=predict_registry(registry,x,bi,cells);expected=apply_policy(pred['primary'],cells,distance,policy);errors=[]
    for i in range(len(chosen)):
        errors.append(abs(float(got[i]['calibrated_laboratory_call_probability'])-expected['calibrated_call_probability'][i]));assert got[i]['status']==str(expected['status'][i]),(i,got[i]['status'],expected['status'][i])
        if got[i]['status']!='ranked':assert got[i]['ranked_score']==''
    assert max(errors)<.002,errors
    write_json(OUT/'scorer_smoke_audit.json',{'passed':True,'actual_CLI_invoked':True,'supported_calibration_inputs':len(chosen),'unknown_cell_abstains':True,'invalid_sequence_abstains':True,'unranked_rows_have_no_ranked_score':True,'cached_vs_fresh_sequence_calibrated_probability_errors':errors,'max_error':max(errors),'predeclared_probability_tolerance':.002,'source_code_sha256':sha(WORK/'score.py'),'scope':'Operational scorer and local numerical reconstruction. These are already observed calibration inputs, not a new biological test.'});print('SCORER OPERATIONAL AUDIT PASSED',max(errors),flush=True)
if __name__=='__main__':main()

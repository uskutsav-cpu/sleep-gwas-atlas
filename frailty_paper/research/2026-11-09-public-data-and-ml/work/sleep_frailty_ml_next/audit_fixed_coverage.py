"""Independent fitted math and exhaustive matched-arm/partition verification."""
import sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
from coverage_inputs import load_inputs
from coverage_validation_common import DEST,FINAL
from matched_coverage_head import MatchedCoverageHead,torch
from safetensors.torch import load_file
from training_utilities import transformations
from activity_effect_network import activity_loss
from audit_first_fitted_network import numpy_loss
from audit_cell_matched_fitted_head import independent_forward

def main():
    start=time.monotonic();assert_previous_sources();settings=json.loads((FINAL/'final_settings_lock.json').read_text());d,pairs,bank,x=load_inputs();checks=[]
    for source in settings['held_studies']:
        train=d['split']!=2
        if source=='none':test=d['split']==2;other=None
        else:
            sid=list(d['study_names']).index(source);groups=np.unique(d['group'][d['study']==sid]);train&=(d['study']!=sid)&~np.isin(d['group'],groups);test=d['study']==sid;other=(d['split']==2)&(d['study']!=sid)&~np.isin(d['group'],groups)
        tr,te=np.flatnonzero(train),np.flatnonzero(test);assert not set(d['group'][tr])&set(d['group'][te])
        for seed in settings['seeds']:
            records=[json.loads((FINAL/'fits'/f'fixed_{source}_{arm}_{seed}_record.json').read_text()) for arm in ['sequence','AVI']]
            for key in ['training_index_sha256','prediction_index_sha256','shared_initial_parameter_sha256','actual_batch_schedule_sha256','completed_epochs','actual_backprop_updates']:assert records[0][key]==records[1][key],(source,seed,key)
            for r in records:
                assert r['settings_lock_sha256']==sha(FINAL/'final_settings_lock.json');assert r['training_index_sha256']==hashlib.sha256(tr.tobytes()).hexdigest();assert r['prediction_index_sha256']==hashlib.sha256(te.tobytes()).hexdigest()
                for name in ['checkpoint','predictions']:
                    assert sha(r[name])==r[name+'_sha256']
                with np.load(r['predictions']) as z:assert np.array_equal(z['observation_index'],te)
                if other is not None:
                    assert sha(r['other_calibration_predictions'])==r['other_calibration_predictions_sha256']
                    with np.load(r['other_calibration_predictions']) as z:assert np.array_equal(z['observation_index'],np.flatnonzero(other))
            checks.append({'source':source,'seed':seed,'matched':True})
    math_results={}
    for arm in ['sequence','AVI']:
        r=json.loads((FINAL/'fits'/f'fixed_none_{arm}_42_record.json').read_text());train=d['split']!=2
        main_ids=np.r_[np.flatnonzero(train&np.isfinite(d['signed_effect']))[:4],np.flatnonzero(train&np.isfinite(d['activity_reference'])&~np.isfinite(d['signed_effect']))[:4]];blocks=[]
        for key in ['local','cross']:blocks.append(pairs[key][np.all(train[pairs[key]],axis=1)][-1])
        ids=np.r_[main_ids,*blocks];b=d['bank_index'][ids];features=np.clip((x[b]-np.array(r['feature_mean'],dtype='float32'))/np.array(r['feature_sd'],dtype='float32'),-12,12).astype('float64')
        rr=bank[b,:,0].astype('float64');aa=bank[b,:,1].astype('float64');c=d['cell'][ids].astype('int64');line=d['line'][ids];cf=np.eye(len(r['cell_medians']))[c]*np.asarray(r['supported_cells'])[c,None];lf=np.eye(len(r['line_medians']))[line]*np.asarray(r['supported_lines'])[line,None];extra=None
        if arm=='AVI':
            from feature_banks import avi_bank
            old=avi_bank();values=np.empty((len(b),3),dtype='float32');mask=b<old.shape[0];values[mask]=old[b[mask]];values[~mask]=np.load(DEST/'coverage_appended_AVI.npy')[b[~mask]-old.shape[0]];t=r['AVI_transform'];values[:,:2]=(np.where(np.isfinite(values[:,:2]),values[:,:2],t['median'])-t['mean'])/t['sd'];extra=np.clip(values,-12,12).astype('float64')
        model=MatchedCoverageHead(r['cell_medians'],r['line_medians'],r['line_parent'],42,arm=='AVI').double().eval();model.load_state_dict(load_file(r['checkpoint']));state={k:v.detach().numpy().copy() for k,v in model.state_dict().items()};nx=np.c_[features,extra] if extra is not None else features;expected=independent_forward(state,nx,rr,aa,c,cf,lf)
        inputs=[torch.from_numpy(v) for v in [features,rr,aa,c,cf,lf]]+[torch.from_numpy(extra) if extra is not None else None];got=model(*inputs);errors={k:float(abs(got[k].detach().numpy()-expected[k]).max()) for k in expected};assert max(errors.values())<1e-7
        med,support,w,centers,quality,context=transformations(d,train,True);assert context==r['activity_and_precision_transformations'];batch={k:d[k][ids].astype('float64') for k in ['effect','signed_effect','activity_reference','activity_alternate']};batch['activity_center']=centers[ids].astype('float64');tb={k:torch.from_numpy(v) for k,v in batch.items()};tw=torch.from_numpy(w[ids].astype('float64'));tq=torch.from_numpy(quality[ids].astype('float64'));loss=activity_loss(got,tb,tw,tq,8,1,1);gold=numpy_loss(expected,batch,w[ids].astype('float64'),quality[ids].astype('float64'),8,1,1);assert abs(float(loss.detach())-gold)<1e-8;loss.backward();gradients=[]
        names=['cell_activity_weight','line_activity_weight','hidden.weight','hidden2.weight','output.weight','cell_bias','line_offset']+(['AVI_added_input_columns'] if arm=='AVI' else [])
        for name in names:
            parameter=dict(model.named_parameters())['hidden.weight' if name=='AVI_added_input_columns' else name]
            if name=='AVI_added_input_columns':local=parameter.grad[:,66:69];flat=int(local.abs().argmax());row,column=divmod(flat,3);index=row*parameter.shape[1]+66+column
            else:index=int(parameter.grad.abs().argmax())
            analytic=float(parameter.grad.reshape(-1)[index]);assert abs(analytic)>1e-12;original=float(parameter.detach().reshape(-1)[index]);epsilon=1e-6
            with torch.no_grad():
                parameter.reshape(-1)[index]=original+epsilon;plus=float(activity_loss(model(*inputs),tb,tw,tq,8,1,1));parameter.reshape(-1)[index]=original-epsilon;minus=float(activity_loss(model(*inputs),tb,tw,tq,8,1,1));parameter.reshape(-1)[index]=original
            finite=(plus-minus)/(2*epsilon);relative=abs(finite-analytic)/max(abs(finite),abs(analytic),1e-12);assert relative<1e-4,(arm,name,relative);gradients.append({'parameter':name,'relative_gradient_error':relative})
        math_results[arm]={'NumPy_forward_maximum_errors':errors,'NumPy_loss_error':abs(float(loss.detach())-gold),'fitted_gradient_checks':gradients,'checkpoint_sha256':r['checkpoint_sha256']}
    write_json(FINAL/'fixed_fits_independent_audit.json',{'passed':True,'all_24_matched_arm_comparisons':checks,'fitted_math':math_results,'settings_lock_sha256':sha(FINAL/'final_settings_lock.json'),
              'all_prediction_and_other_calibration_ids_reconstructed':True,'code_sha256':sha(Path(__file__)),'scope':'Numerical, matching, partition and source-omission correctness. No assertion of useful prediction or causal disease biology.','elapsed_seconds':time.monotonic()-start});assert_previous_sources();print('FIXED COVERAGE FITS INDEPENDENT AUDIT PASSED',flush=True)

if __name__=='__main__':main()

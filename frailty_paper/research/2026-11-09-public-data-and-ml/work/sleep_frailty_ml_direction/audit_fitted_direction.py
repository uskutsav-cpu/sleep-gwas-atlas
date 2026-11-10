"""Audit fitted gradients on actual newly recovered signed examples."""
from direction_common import *
from measurement_inputs import load_data, normalized_inputs
from folded_effect_network import FoldedEffectHead, torch
from audit_fitted_measurement import numpy_forward, numpy_loss
from activity_effect_network import activity_loss
from safetensors.torch import load_file


def main():
    assert_versions();d,pairs,bank,x,extra=load_data();new=load_signed_sidecar(d);train=d['split']==0
    main=[]
    for c in np.unique(d['line'][new&train]):main.extend(np.flatnonzero(new&train&(d['line']==c))[:2])
    kosicki=list(d['study_names']).index('Kosicki2025');main.extend(np.flatnonzero(train&(d['study']==kosicki)&np.isfinite(d['signed_effect']))[:4])
    pair_blocks=[p[np.all(train[p],axis=1)][0] for p in [pairs['local'],pairs['cross']]]
    ids=np.r_[main,*pair_blocks];nmain=len(main);assert new[main].sum()==8
    stem=DEST/'development/source_signed_none_42';record=json.loads(Path(str(stem)+'_record.json').read_text())
    checkpoint=Path(str(stem)+'.safetensors');assert sha(checkpoint)==record['checkpoint_sha256']
    b=d['bank_index'][ids];mean,difference,physical=[a.astype('float64') for a in normalized_inputs(extra,x,b,record['normalization'])]
    raw=bank[b].astype('float64');cell=d['cell'][ids].astype('int64');line=d['line'][ids]
    cf=np.eye(len(record['cell_medians']))[cell]*np.array(record['supported_cells'])[cell,None]
    lf=np.eye(len(record['line_medians']))[line]*np.array(record['supported_lines'])[line,None]
    values=[mean,difference,physical,raw[:,:,0],raw[:,:,1],cell,cf,lf];tensors=[torch.from_numpy(a.copy()) for a in values]
    from training_utilities import transformations
    _,_,weights,centers,quality,_=transformations(d,train,True)
    batch={k:d[k][ids].astype('float64') for k in ['effect','signed_effect','activity_reference','activity_alternate']}
    batch['activity_center']=centers[ids].astype('float64');weights=weights[ids].astype('float64');quality=quality[ids].astype('float64')
    tb={k:torch.from_numpy(v) for k,v in batch.items()};tw,tq=map(torch.from_numpy,[weights,quality])
    model=FoldedEffectHead(record['cell_medians'],record['line_medians'],record['line_parent']).double().eval()
    model.load_state_dict(load_file(str(checkpoint)),strict=True);state={k:v.detach().numpy().copy() for k,v in model.state_dict().items()}
    expected=numpy_forward(state,*values);got=model(*tensors)
    errors={k:float(np.max(abs(got[k].detach().numpy()-expected[k]))) for k in expected};assert max(errors.values())<1e-8
    loss=activity_loss(got,tb,tw,tq,nmain,1,1);independent=numpy_loss(expected,batch,weights,quality,nmain,1,1)
    assert abs(float(loss.detach())-independent)<1e-8
    loss.backward();checks=[]
    for name in ['delta_hidden.weight','delta_hidden2.weight','delta_output.weight','noise_hidden.weight','noise_output.weight']:
        parameter=dict(model.named_parameters())[name];index=int(parameter.grad.abs().argmax());analytic=float(parameter.grad.reshape(-1)[index]);assert abs(analytic)>1e-10
        original=float(parameter.detach().reshape(-1)[index]);epsilon=2e-6;objectives=[]
        for sign in [1,-1]:
            altered={k:v.copy() for k,v in state.items()};altered[name].reshape(-1)[index]=original+sign*epsilon
            objectives.append(numpy_loss(numpy_forward(altered,*values),batch,weights,quality,nmain,1,1))
        finite=(objectives[0]-objectives[1])/(2*epsilon);relative=abs(finite-analytic)/max(abs(finite),abs(analytic),1e-12)
        assert relative<1e-4,(name,analytic,finite,relative);checks.append({'parameter':name,'relative_gradient_error':relative})
    control=load_file(str(MEASUREMENT/'development/observation_only_none_42.safetensors'))
    differences={k:float(np.max(abs(state[k]-v.numpy()))) for k,v in control.items()}
    assert differences['delta_hidden.weight']>0 and differences['delta_output.weight']>0
    assert sha(checkpoint)==record['checkpoint_sha256'];assert_versions()
    write_json(DEST/'fitted_direction_independent_audit.json',{'passed':True,
       'actual_new_signed_rows_checked':8,'new_signed_training_rows_used_in_actual_fit':record['new_original_signed_training_targets'],
       'NumPy_forward_errors':errors,'NumPy_loss_error':abs(float(loss.detach())-independent),
       'fitted_gradients':checks,'maximum_parameter_disagreement_vs_control':differences,
       'checkpoint_sha256':record['checkpoint_sha256'],'saved_weights_changed_by_audit':False,
       'gradient_scope':'Signed and noise branches only; predicted activity is intentionally detached as a noise conditioner.',
       'new_independent_prediction_success_established':False,'reserved_study_outcomes_read':False})
    print('FITTED SOURCE-DIRECTION AUDIT PASS',flush=True)


if __name__=='__main__':main()

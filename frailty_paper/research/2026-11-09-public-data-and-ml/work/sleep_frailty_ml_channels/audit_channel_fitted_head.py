"""Independent NumPy/loss replay and fitted new-channel gradient check."""
from channel_common import *
from coverage_inputs import load_inputs
from training_utilities import transformations
from cell_matched_effect_network import CellMatchedEffectHead,torch
from safetensors.torch import load_file
from activity_effect_network import activity_loss
from audit_first_fitted_network import numpy_loss
from audit_cell_matched_fitted_head import independent_forward

def main():
    assert_frozen_previous();d,pairs,bank,x=load_inputs();extra=np.load(DEST/'projected_reader_channels.npy',mmap_mode='r');train=d['split']==0;main_ids=np.r_[np.flatnonzero(train&np.isfinite(d['signed_effect']))[:4],np.flatnonzero(train&np.isfinite(d['activity_reference'])&~np.isfinite(d['signed_effect']))[:4]];blocks=[]
    for label in ['local','cross']:blocks.append(pairs[label][np.all(train[pairs[label]],axis=1)][-1])
    ids=np.r_[main_ids,*blocks];b=d['bank_index'][ids];results={}
    for recipe in ['scalar_mask','channel_preserving']:
        record=DEST/'development'/f'{recipe}_none_42_record.json';r=json.loads(record.read_text());checkpoint=DEST/'development'/f'{recipe}_none_42.safetensors';assert sha(checkpoint)==r['checkpoint_sha256']
        base=np.clip((x[b]-np.asarray(r['scalar_feature_mean'],dtype='float32'))/np.asarray(r['scalar_feature_sd'],dtype='float32'),-12,12);added=np.clip((extra[b].reshape(-1,192).astype('float32')-np.asarray(r['channel_feature_mean'],dtype='float32'))/np.asarray(r['channel_feature_sd'],dtype='float32'),-12,12)
        if recipe=='scalar_mask':added[:]=0
        features=np.c_[base,added].astype('float64');rr=bank[b,:,0].astype('float64');aa=bank[b,:,1].astype('float64');c=d['cell'][ids].astype('int64');line=d['line'][ids];cf=np.eye(len(r['cell_medians']))[c]*np.asarray(r['supported_cells'])[c,None];lf=np.eye(len(r['line_medians']))[line]*np.asarray(r['supported_lines'])[line,None]
        model=CellMatchedEffectHead(r['cell_medians'],r['line_medians'],r['line_parent'],258,42).double().eval();model.load_state_dict(load_file(checkpoint));state={k:v.detach().numpy().copy() for k,v in model.state_dict().items()};expected=independent_forward(state,features,rr,aa,c,cf,lf);inputs=[torch.from_numpy(v) for v in [features,rr,aa,c,cf,lf]];got=model(*inputs);errors={k:float(abs(got[k].detach().numpy()-expected[k]).max()) for k in expected};assert max(errors.values())<1e-7
        med,support,w,centers,quality,context=transformations(d,train,True);assert context==r['context_transformations'];batch={k:d[k][ids].astype('float64') for k in ['effect','signed_effect','activity_reference','activity_alternate']};batch['activity_center']=centers[ids].astype('float64');tb={k:torch.from_numpy(v) for k,v in batch.items()};tw=torch.from_numpy(w[ids].astype('float64'));tq=torch.from_numpy(quality[ids].astype('float64'));loss=activity_loss(got,tb,tw,tq,8,1,1);gold=numpy_loss(expected,batch,w[ids].astype('float64'),quality[ids].astype('float64'),8,1,1);assert abs(float(loss.detach())-gold)<1e-8;loss.backward();checks=[]
        hidden=dict(model.named_parameters())['hidden.weight'];channel_gradient=hidden.grad[:,66:258];assert torch.isfinite(channel_gradient).all()
        if recipe=='scalar_mask':assert torch.count_nonzero(channel_gradient)==0
        else:assert torch.count_nonzero(channel_gradient)>0
        for name in ['hidden.weight','hidden2.weight','output.weight','cell_bias','line_offset']+(['added_channel_columns'] if recipe=='channel_preserving' else []):
            parameter=hidden if name=='added_channel_columns' else dict(model.named_parameters())[name]
            if name=='added_channel_columns':flat=int(channel_gradient.abs().argmax());row,column=divmod(flat,192);index=row*hidden.shape[1]+66+column
            else:index=int(parameter.grad.abs().argmax())
            analytic=float(parameter.grad.reshape(-1)[index]);assert abs(analytic)>1e-12;original=float(parameter.detach().reshape(-1)[index]);epsilon=1e-6
            with torch.no_grad():
                parameter.reshape(-1)[index]=original+epsilon;plus=float(activity_loss(model(*inputs),tb,tw,tq,8,1,1));parameter.reshape(-1)[index]=original-epsilon;minus=float(activity_loss(model(*inputs),tb,tw,tq,8,1,1));parameter.reshape(-1)[index]=original
            finite=(plus-minus)/(2*epsilon);relative=abs(finite-analytic)/max(abs(finite),abs(analytic),1e-12);assert relative<1e-4,(recipe,name,relative);checks.append({'parameter':name,'relative_gradient_error':relative})
        results[recipe]={'forward_maximum_errors':errors,'independent_loss_error':abs(float(loss.detach())-gold),'gradient_checks':checks,'additional_channel_nonzero_gradient_entries':int(torch.count_nonzero(channel_gradient)),'checkpoint_sha256':sha(checkpoint)}
    write_json(DEST/'fitted_channel_independent_audit.json',{'passed':True,'recipes':results,'prior_frozen_files_unchanged':True,'scientific_usefulness_or_new_validation_established':False,'code_sha256':sha(Path(__file__))});assert_frozen_previous();print('CHANNEL FITTED MATH AND NEW GRADIENTS VERIFIED',flush=True)

if __name__=='__main__':main()

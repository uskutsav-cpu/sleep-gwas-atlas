"""Independent NumPy forward/loss and finite-difference gradients for the new head."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
from training_utilities import transformations
from audit_first_fitted_network import numpy_loss
from activity_effect_network import activity_loss
from cell_matched_effect_network import CellMatchedEffectHead,torch
from safetensors.torch import load_file

DEST=ROOT/'outputs/sleep_frailty_ml_next'


def independent_forward(s,x,rr,aa,c,cf,lf):
    silu=lambda v:v/(1+np.exp(np.clip(-v,-600,600)))
    weights=s['cell_activity_weight'][c]+lf@s['line_activity_weight'];bias=s['cell_activity_bias'][c]+lf@s['line_activity_bias'];reference=(rr*weights).sum(1)+bias;alternate=(aa*weights).sum(1)+bias
    inputs=np.c_[x,abs(alternate-reference)*10.,(reference+alternate)*.5,cf,lf];h=silu(inputs@s['hidden.weight'].T+s['hidden.bias']);h=silu(h@s['hidden2.weight'].T+s['hidden2.bias'])
    raw=s['cell_bias'][c]+lf@s['line_offset']+h@s['output.weight'][0]+s['output.bias'][0]
    return {'magnitude':np.logaddexp(0,raw),'signed_effect':alternate-reference,'activity_reference':reference,'activity_alternate':alternate}


def main():
    assert_previous_sources();r=json.loads((DEST/'development/three_readers_none_42_record.json').read_text());checkpoint=DEST/'development/three_readers_none_42.safetensors';assert sha(checkpoint)==r['checkpoint_sha256']
    with np.load(DATA/'training_arrays.npz') as z:d={k:z[k] for k in z.files}
    with np.load(DATA/'measured_pairs.npz') as z:pairs={k:z[k] for k in z.files}
    train=d['split']==0;main_ids=np.r_[np.flatnonzero(train&np.isfinite(d['signed_effect']))[:4],np.flatnonzero(train&np.isfinite(d['activity_reference'])&~np.isfinite(d['signed_effect']))[:4]];blocks=[]
    for key in ['local','cross']:
        q=pairs[key][np.all(train[pairs[key]],axis=1)];blocks.append(q[-1])
    ids=np.r_[main_ids,*blocks];assert len(ids)==12;b=d['bank_index'][ids];bank=np.load(DEST/'cell_matched_native_development_features.npy',mmap_mode='r');pp=np.load(DEST/'native_development_physical_metadata.npy')[b]
    x=np.c_[bank[b].reshape(12,48),pp[:,0]/230.,pp[:,1]/np.maximum(pp[:,0]-1,1),np.eye(16)[pp[:,2]*4+pp[:,3]]].astype('float32');x=np.clip((x-np.array(r['feature_mean'],dtype='float32'))/np.array(r['feature_sd'],dtype='float32'),-12,12).astype('float64')
    rr=bank[b,:,0].astype('float64');aa=bank[b,:,1].astype('float64');c=d['cell'][ids].astype('int64');line=d['line'][ids];cf=np.eye(len(r['cell_medians']))[c]*np.array(r['supported_cells'])[c,None];lf=np.eye(len(r['line_medians']))[line]*np.array(r['supported_lines'])[line,None]
    model=CellMatchedEffectHead(r['cell_medians'],r['line_medians'],r['line_parent']).double().eval();model.load_state_dict(load_file(checkpoint));s={k:v.detach().numpy().copy() for k,v in model.state_dict().items()};expected=independent_forward(s,x,rr,aa,c,cf,lf)
    inputs=[torch.from_numpy(v) for v in [x,rr,aa,c,cf,lf]];got=model(*inputs);errors={k:float(abs(got[k].detach().numpy()-expected[k]).max()) for k in expected};assert max(errors.values())<1e-7,errors
    med,support,w,centers,quality,context=transformations(d,train,True);assert context==r['context_transformations'];batch={k:d[k][ids].astype('float64') for k in ['effect','signed_effect','activity_reference','activity_alternate']};batch['activity_center']=centers[ids].astype('float64');tb={k:torch.from_numpy(v) for k,v in batch.items()};tw=torch.from_numpy(w[ids].astype('float64'));tq=torch.from_numpy(quality[ids].astype('float64'))
    loss=activity_loss(got,tb,tw,tq,8,1,1);gold=numpy_loss(expected,batch,w[ids].astype('float64'),quality[ids].astype('float64'),8,1,1);assert abs(float(loss.detach())-gold)<1e-8;loss.backward();checks=[]
    for name in ['cell_activity_weight','line_activity_weight','hidden.weight','hidden2.weight','output.weight','cell_bias','line_offset']:
        p=dict(model.named_parameters())[name];assert p.grad is not None and torch.isfinite(p.grad).all();index=int(p.grad.abs().argmax());analytic=float(p.grad.reshape(-1)[index]);assert abs(analytic)>1e-12
        original=float(p.detach().reshape(-1)[index]);epsilon=1e-6
        with torch.no_grad():
            p.reshape(-1)[index]=original+epsilon;plus=float(activity_loss(model(*inputs),tb,tw,tq,8,1,1));p.reshape(-1)[index]=original-epsilon;minus=float(activity_loss(model(*inputs),tb,tw,tq,8,1,1));p.reshape(-1)[index]=original
        finite=(plus-minus)/(2*epsilon);relative=abs(finite-analytic)/max(abs(finite),abs(analytic),1e-12);assert relative<1e-4,(name,relative);checks.append({'parameter':name,'relative_gradient_error':relative})
    write_json(DEST/'fitted_head_independent_audit.json',{'passed':True,'NumPy_forward_maximum_errors':errors,'independent_NumPy_loss_error':abs(float(loss.detach())-gold),'fitted_gradient_checks':checks,'checkpoint_sha256':sha(checkpoint),'new_test_cases_or_outcomes_used':False,'scope':'Computation and fitted gradients only. No claim of a useful molecular or disease predictor.'});assert_previous_sources();print('CELL-MATCHED FITTED HEAD INDEPENDENTLY VERIFIED',flush=True)


if __name__=='__main__':main()

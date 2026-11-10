"""Actual native-DNA, matched-state, finite-gradient and update audit."""
import time
from native_common import *
from native_data import load_data,PackedDNA,normalization
from native_network import NativeAdaptationModel,torch,author_reader,local_pool,summarize
from activity_effect_network import activity_loss
from training_utilities import transformations
from extract_cell_matched_development_features import summarize as numpy_summary


def main():
    assert_versions_native();plan=get_plan();start=time.monotonic();d,pairs,bank,x,extra=load_data();first,neural=add_all_documented_signed_targets(d)
    packed=PackedDNA(d,bank);train=d['split']==0;tr=np.flatnonzero(train);med,support,w,centers,quality,ctx=transformations(d,train,True)
    nl=len(d['cell_line_names']);parent=np.zeros(nl,dtype='int64');lm=np.zeros(nl,dtype='float32');ls=np.zeros(nl,dtype=bool)
    for line in range(nl):
        cell=np.unique(d['cell'][d['line']==line]);parent[line]=cell[0] if len(cell) else 0
        good=train&(d['line']==line);ls[line]=good.any();lm[line]=np.median(d['effect'][good]) if good.any() else med[parent[line]]
    norm=normalization(bank,x,extra,np.unique(d['bank_index'][tr]));chosen=[]
    for length in np.unique(packed.length[np.unique(d['bank_index'][tr])]):
        chosen.extend(tr[packed.length[d['bank_index'][tr]]==length][:2])
    chosen.extend(np.flatnonzero(train&neural)[:4]);chosen.extend(np.flatnonzero(train&np.isfinite(d['signed_effect'])&~first&~neural)[:4])
    ids=np.array(chosen,dtype='int64');assert len(ids)>=12 and neural[ids].any()
    torch.manual_seed(42);frozen=NativeAdaptationModel(med,lm,parent,norm,'frozen_reader').double().eval()
    torch.manual_seed(42);adapt=NativeAdaptationModel(med,lm,parent,norm,'adapt_last_stage').double().eval()
    assert all(torch.equal(v,adapt.state_dict()[k]) for k,v in frozen.state_dict().items())
    a=frozen(ids,d,packed,bank,x,extra,support,ls);b=adapt(ids,d,packed,bank,x,extra,support,ls)
    assert all(torch.equal(a[k],b[k]) for k in a)
    reverse=adapt(ids,d,packed,bank,x,extra,support,ls,True)
    assert torch.allclose(b['magnitude'],reverse['magnitude'],atol=1e-9,rtol=0)
    assert torch.allclose(b['signed_effect'],-reverse['signed_effect'],atol=1e-9,rtol=0)
    # Independent NumPy summary replay from actual native activations.
    which=np.array([int(d['bank_index'][ids[0]])]);tokens=packed.batch(which,'cpu',torch.float64)
    with torch.no_grad():
        acts=adapt.reader.main(adapt.reader.stem(tokens));g=acts.mean(-1);position=torch.tensor(np.r_[packed.position[which],packed.position[which]],dtype=torch.int64)
        local=local_pool(acts,position);p=adapt.reader.head(adapt.reader.mapper(acts).mean(-1)).squeeze(-1)
        summaries,blocks=summarize(g,local,p,1);expected=numpy_summary(g.numpy(),local.numpy(),p.numpy(),1)
        summary_error=float(abs(expected-summaries.numpy()).max())
        # The pinned independent summary function returns float32; compare to
        # its actual storage precision rather than claiming double precision.
        assert summary_error<2*np.finfo(np.float32).eps*max(1.,float(abs(expected).max())),summary_error
    batch={k:torch.from_numpy(d[k][ids]).double() for k in ['effect','signed_effect','activity_reference','activity_alternate']}
    batch['activity_center']=torch.from_numpy(centers[ids]).double();tw=torch.from_numpy(w[ids]).double();tq=torch.from_numpy(quality[ids]).double()
    def objective(model):return activity_loss(model(ids,d,packed,bank,x,extra,support,ls),batch,tw,tq,len(ids),0,0)+model.anchor_loss()
    loss=objective(adapt);loss.backward();checks=[]
    for name in ['reader.main.3.0.fn.block.0.weight','reader.main.3.1.block.0.weight']:
        parameter=dict(adapt.named_parameters())[name];assert parameter.grad is not None and torch.isfinite(parameter.grad).all()
        index=int(parameter.grad.abs().argmax());analytic=float(parameter.grad.reshape(-1)[index]);assert abs(analytic)>1e-9
        original=float(parameter.detach().reshape(-1)[index]);epsilon=1e-5
        with torch.no_grad():
            parameter.reshape(-1)[index]=original+epsilon;plus=float(objective(adapt))
            parameter.reshape(-1)[index]=original-epsilon;minus=float(objective(adapt));parameter.reshape(-1)[index]=original
        finite=(plus-minus)/(2*epsilon);relative=abs(finite-analytic)/max(abs(finite),abs(analytic),1e-12)
        assert relative<5e-4,(name,finite,analytic,relative);checks.append({'parameter':name,'relative_gradient_error':relative})
    device='mps' if torch.backends.mps.is_available() else 'cpu'
    torch.manual_seed(42);model=NativeAdaptationModel(med,lm,parent,norm,'adapt_last_stage').to(device).eval()
    initial={k:v.detach().cpu().clone() for k,v in model.reader.state_dict().items()}
    with torch.no_grad():gpu=model(ids,d,packed,bank,x,extra,support,ls)
    errors={k:float(abs(gpu[k].cpu().numpy()-b[k].detach().numpy()).max()) for k in b};assert max(errors.values())<.003,errors
    optimizer=torch.optim.AdamW([{'params':model.head.parameters(),'lr':3e-4},
        {'params':[p for p in model.reader.parameters() if p.requires_grad],'lr':1e-5}],weight_decay=.001)
    tb={k:v.to(device=device,dtype=torch.float32) for k,v in batch.items()};weights=tw.to(device=device,dtype=torch.float32);q=tq.to(device=device,dtype=torch.float32)
    losses=[];step_start=time.monotonic()
    for step in range(3):
        optimizer.zero_grad();prediction=model(ids,d,packed,bank,x,extra,support,ls)
        value=activity_loss(prediction,tb,weights,q,len(ids),0,0)+model.anchor_loss();assert torch.isfinite(value)
        value.backward();assert all(p.grad is None or torch.isfinite(p.grad).all() for p in model.parameters())
        torch.nn.utils.clip_grad_norm_([p for p in model.parameters() if p.requires_grad],1.);optimizer.step();losses.append(float(value.detach().cpu()))
    if device=='mps':torch.mps.synchronize()
    changes={k:int(torch.count_nonzero(v.detach().cpu()!=initial[k])) for k,v in model.reader.state_dict().items()}
    assert sum(v for k,v in changes.items() if k.startswith('main.3.'))>0
    assert sum(v for k,v in changes.items() if not k.startswith('main.3.'))==0
    assert sum(v for k,v in changes.items() if 'running_' in k or 'num_batches' in k)==0
    after=model(ids,d,packed,bank,x,extra,support,ls);swapped=model(ids,d,packed,bank,x,extra,support,ls,True)
    assert torch.allclose(after['magnitude'],swapped['magnitude'],atol=1e-6,rtol=0)
    assert torch.allclose(after['signed_effect'],-swapped['signed_effect'],atol=1e-6,rtol=0)
    write_json(DEST/'native_runtime_and_gradient_audit.json',{'passed':True,'device':device,'actual_training_observation_indices':ids.tolist(),
       'actual_new_neuronal_signed_labels_in_audit':int(neural[ids].sum()),'packed_physical_DNA':packed.receipt,
       'same_initial_parameter_values_and_outputs_between_arms':True,'magnitude_swap_invariance_and_signed_antisymmetry_before_and_after_updates':True,
       'independent_NumPy_summary_max_error':summary_error,'independent_summary_reference_dtype':'float32',
       'native_last_stage_finite_gradient_checks':checks,
       'CPU_double_to_runtime_float32_output_errors':errors,'real_backprop_updates':3,'engineering_losses':losses,
       'changed_reader_scalar_counts':changes,'only_last_stage_reader_parameters_changed':True,'BN_running_statistics_changed':False,
       'engineering_updates_not_used_to_initialize_scientific_fits':True,'small_batch_three_step_seconds':time.monotonic()-step_start,
       'new_outside_test_read':False,'scientific_predictive_improvement_established':False,'elapsed_seconds':time.monotonic()-start})
    assert_versions_native();print('NATIVE DNA READER AND ACTUAL LAST-STAGE GRADIENT AUDIT PASS',flush=True)


if __name__=='__main__':main()

"""Replay a saved native-adaptation fit without changing any saved weights."""
import argparse
from bounded_common import *
from native_data import load_data,PackedDNA
from bounded_network import NativeAdaptationModel,torch
from stable_loss import activity_loss
from training_utilities import transformations
from safetensors.torch import load_file
from audit_first_fitted_network import numpy_loss


def numpy_head(state,args):
    x,rr,aa,cell,cf,lf=args
    def silu(a):return a/(1+np.exp(np.clip(-a,-600,600)))
    weight=state['cell_activity_weight'][cell]+lf@state['line_activity_weight'];bias=state['cell_activity_bias'][cell]+lf@state['line_activity_bias']
    reference=(rr*weight).sum(1)+bias;alternate=(aa*weight).sum(1)+bias;signed=alternate-reference
    z=np.c_[x,abs(signed)*10,(reference+alternate)*.5,cf,lf]
    hidden=silu(z@state['hidden.weight'].T+state['hidden.bias']);hidden=silu(hidden@state['hidden2.weight'].T+state['hidden2.bias'])
    magnitude=np.logaddexp(0,state['cell_bias'][cell]+lf@state['line_offset']+hidden@state['output.weight'][0]+state['output.bias'][0])
    return {'magnitude':magnitude,'signed_effect':signed,'activity_reference':reference,'activity_alternate':alternate}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--recipe',default='adapt_last_stage');parser.add_argument('--source',default='none');parser.add_argument('--epoch',type=int,default=1);args=parser.parse_args()
    assert_versions_stable();stem=f'{args.recipe}_{args.source}_42';directory=DEST/'development'
    meta=json.loads((directory/(stem+'_started.json')).read_text());receipt=json.loads((directory/f'{stem}_epoch{args.epoch}_receipt.json').read_text())
    path=directory/f'{stem}_epoch{args.epoch}.safetensors';assert sha(path)==receipt['checkpoint_sha256']
    d,pairs,bank,x,extra=load_data();first,neural=add_all_documented_signed_targets(d);packed=PackedDNA(d,bank);train=d['split']==0
    if args.source!='none':
        sid=list(d['study_names']).index(args.source);groups=np.unique(d['group'][d['study']==sid]);train&=(d['study']!=sid)&~np.isin(d['group'],groups)
    assert hashlib.sha256(np.flatnonzero(train).tobytes()).hexdigest()==meta['training_index_sha256']
    main=np.r_[np.flatnonzero(train&neural)[:4],np.flatnonzero(train&first)[:4],np.flatnonzero(train&np.isfinite(d['signed_effect'])&~neural&~first)[:4]]
    blocks=[p[np.all(train[p],axis=1)][0] for p in [pairs['local'],pairs['cross']]];ids=np.r_[main,*blocks];nmain=len(main)
    model=NativeAdaptationModel(meta['cell_medians'],meta['line_medians'],meta['line_parent'],meta['normalization'],args.recipe).double().eval()
    state=load_file(str(path));assert set(state)==set(receipt['parameter_names']);missing,unexpected=model.load_state_dict(state,strict=False);assert not unexpected
    # Overrides contain learned parameters, while unchanged BN buffers and the
    # fixed prefix are restored from the pinned original reader.
    parameter_names=set(dict(model.named_parameters()))
    assert all(not n.startswith('head.') and not (args.recipe=='adapt_last_stage' and n.startswith('reader.main.3.')) for n in missing if n in parameter_names)
    supported=np.array(meta['supported_cells']);lines=np.array(meta['supported_lines']);captured=[]
    hook=model.head.register_forward_pre_hook(lambda module,inputs:captured.append([v.detach().numpy().copy() for v in inputs]))
    got=model(ids,d,packed,bank,x,extra,supported,lines);hook.remove()
    expected=numpy_head({k:v.detach().numpy().copy() for k,v in model.head.state_dict().items()},captured[0])
    errors={k:float(abs(expected[k]-got[k].detach().numpy()).max()) for k in expected};assert max(errors.values())<1e-8
    _,_,w,centers,q,_=transformations(d,train,True);batch={k:d[k][ids].astype('float64') for k in ['effect','signed_effect','activity_reference','activity_alternate']};batch['activity_center']=centers[ids].astype('float64')
    weights=w[ids].astype('float64');quality=q[ids].astype('float64');tb={k:torch.from_numpy(v) for k,v in batch.items()};tw,tq=map(torch.from_numpy,[weights,quality])
    loss=activity_loss(got,tb,tw,tq,nmain,1,1);independent=numpy_loss(expected,batch,weights,quality,nmain,1,1)
    assert abs(float(loss.detach())-independent)<1e-8
    def objective():return activity_loss(model(ids,d,packed,bank,x,extra,supported,lines),tb,tw,tq,nmain,1,1)+model.anchor_loss()
    checks=[]
    if args.recipe=='adapt_last_stage':
        full=objective();full.backward()
        for name in ['reader.main.3.0.fn.block.0.weight','reader.main.3.1.block.0.weight']:
            p=dict(model.named_parameters())[name];assert p.grad is not None and torch.isfinite(p.grad).all();index=int(p.grad.abs().argmax());analytic=float(p.grad.reshape(-1)[index]);assert abs(analytic)>1e-10
            original=float(p.detach().reshape(-1)[index]);eps=1e-6
            with torch.no_grad():
                p.reshape(-1)[index]=original+eps;plus=float(objective());p.reshape(-1)[index]=original-eps;minus=float(objective());p.reshape(-1)[index]=original
            finite=(plus-minus)/(2*eps);relative=abs(finite-analytic)/max(abs(finite),abs(analytic),1e-12)
            assert relative<5e-4,(name,finite,analytic,relative);checks.append({'parameter':name,'relative_gradient_error':relative})
    forward=model(ids,d,packed,bank,x,extra,supported,lines);reverse=model(ids,d,packed,bank,x,extra,supported,lines,True)
    assert torch.allclose(forward['magnitude'],reverse['magnitude'],atol=1e-9,rtol=0)
    assert torch.allclose(forward['signed_effect'],-reverse['signed_effect'],atol=1e-9,rtol=0)
    prediction_path=directory/f'{stem}_epoch{args.epoch}_predictions.npz'
    assert sha(prediction_path)==receipt['prediction_sha256']
    with np.load(prediction_path,allow_pickle=False) as saved:
        te=saved['observation_index'];assert hashlib.sha256(te.tobytes()).hexdigest()==meta['prediction_index_sha256']
        lengths=packed.length[d['bank_index'][te]]
        positions=np.concatenate([np.flatnonzero(lengths==length)[:2] for length in np.unique(lengths)])
        replay_ids=te[positions]
        with torch.no_grad():replayed=model(replay_ids,d,packed,bank,x,extra,supported,lines)
        replay_errors={k:float(abs(replayed[k].cpu().numpy()-saved[k][positions]).max()) for k in replayed}
    # The saved inference used MPS float32; the independent replay is CPU float64.
    # Compare a stated absolute precision tolerance and report the actual error.
    assert max(replay_errors.values())<1e-4,replay_errors
    assert sha(path)==receipt['checkpoint_sha256'];assert_versions_stable()
    write_json(DEST/f'{stem}_epoch{args.epoch}_independent_audit.json',{'passed':True,'checkpoint_sha256':receipt['checkpoint_sha256'],
       'independent_NumPy_head_errors':errors,'independent_NumPy_data_loss_error':abs(float(loss.detach())-independent),
       'actual_fitted_core_finite_gradient_checks':checks,'fitted_allele_swap_size_invariance_and_direction_antisymmetry':True,
       'actual_training_rows_checked':len(ids),'saved_parameters_changed_by_audit':False,
       'saved_validation_prediction_rows_replayed':replay_ids.tolist(),
       'saved_validation_prediction_lengths_checked':np.unique(lengths).tolist(),
       'CPU_double_vs_saved_MPS_float32_max_errors':replay_errors,'replay_absolute_tolerance':1e-4,
       'auditor_code_sha256':sha(Path(__file__)),
       'reserved_new_study_read':False,'scientific_predictive_success_established':False})
    print('FITTED NATIVE MODEL INDEPENDENT AUDIT PASS',stem,args.epoch,flush=True)


if __name__=='__main__':main()

"""Independent NumPy replay and fitted gradients for explicit culture contexts."""
from v9_common import *
from line_aware_network import LineAwareEffectNetwork,torch
from feature_banks import global_bank,local_bank
from training_utilities import transformations
from activity_effect_network import activity_loss
from safetensors.torch import load_file


def numpy_forward(s,g,l,c,cf,lf):
    def bn(x,p):return (x-s[p+'.running_mean'])/np.sqrt(s[p+'.running_var']+1e-5)*s[p+'.weight']+s[p+'.bias']
    def silu(x):return x/(1+np.exp(np.clip(-x,-600,600)))
    def embed(x):
        x=bn(x,'reader.mapper.block.0')@s['reader.mapper.block.1.weight'][:,:,0].T+s['reader.mapper.block.1.bias']
        return silu(bn(x@s['reader.head.0.weight'].T+s['reader.head.0.bias'],'reader.head.1'))
    a,b=embed(g[:,0]),embed(g[:,1]);la,lb=embed(l[:,0]),embed(l[:,1]);w=s['cell_activity.weight'][c]+lf@s['line_activity.weight'];bias=s['cell_activity.bias'][c]+lf@s['line_activity.bias']
    ar=(a*w).sum(1)+bias;aa=(b*w).sum(1)+bias;mean=(a+b)*.5;centered=mean-mean.mean(1,keepdims=True);normalized=centered/np.sqrt((centered**2).mean(1,keepdims=True)+1e-5)
    x=np.c_[abs(b-a)*10.,normalized,abs(aa-ar)*10.,(ar+aa)*.5,np.log1p(np.sqrt((mean*mean).mean(1)+1e-12)),abs(lb-la)*10.,cf,lf]
    h=silu(x@s['mag_hidden.weight'].T+s['mag_hidden.bias']);raw=s['mag_cell_bias'][c]+lf@s['mag_line_offset']+h@s['mag_out.weight'][0]+s['mag_out.bias'][0]
    return {'magnitude':np.logaddexp(0,raw),'signed_effect':aa-ar,'activity_reference':ar,'activity_alternate':aa}


def numpy_loss(p,b,w,q,nmain,nlocal,ncross):
    result=np.mean(w[:nmain]*abs(p['magnitude'][:nmain]-b['effect'][:nmain]))
    signed=np.isfinite(b['signed_effect']);signed[nmain:]=False
    if signed.any():result+=.25*np.mean(w[signed]*q[signed]*abs(p['signed_effect'][signed]-b['signed_effect'][signed]))
    active=np.isfinite(b['activity_reference'])&np.isfinite(b['activity_alternate']);active[nmain:]=False
    def huber(e):a=abs(e);return np.where(a<1,.5*a*a,a-.5)
    if active.any():
        ar=b['activity_reference'][active]-b['activity_center'][active];aa=b['activity_alternate'][active]-b['activity_center'][active]
        result+=.05*np.mean(w[active]*(huber(p['activity_reference'][active]-ar)+huber(p['activity_alternate'][active]-aa))*.5)
        result+=.05*np.mean(w[active]*q[active]*huber(p['signed_effect'][active]-(aa-ar)))
    for count,offset in [(nlocal,nmain),(ncross,nmain+2*nlocal)]:
        pp=p['magnitude'][offset:offset+2*count].reshape(count,2);yy=b['effect'][offset:offset+2*count].reshape(count,2);qq=q[offset:offset+2*count].reshape(count,2)
        result+=.2*np.mean(np.sqrt(qq.prod(1))*abs((pp[:,0]-pp[:,1])-(yy[:,0]-yy[:,1])))
    return float(result)


def main():
    r=json.loads((OUT/'frozen_none_42_record.json').read_text());assert sha(r['checkpoint'])==r['checkpoint_sha256']
    with np.load(DATA/'training_arrays.npz') as z:d={k:z[k] for k in z.files}
    with np.load(DATA/'measured_pairs.npz') as z:pairs={k:z[k] for k in z.files}
    train=d['split']==0;ret=list(d['study_names']).index('Retallick2026');mc=list(d['study_names']).index('McAfee2023')
    main=np.r_[np.flatnonzero(train&np.isfinite(d['signed_effect']))[:4],np.flatnonzero(train&(d['study']==ret)&np.isfinite(d['activity_reference']))[:4]]
    blocks=[]
    for key in ['local','cross']:
        pp=pairs[key];good=np.all(train[pp],axis=1)&(d['study'][pp[:,0]]==mc);assert good.any();blocks.append(pp[good][0])
    ids=np.r_[main,*blocks];assert len(ids)==12
    g=global_bank()[d['bank_index'][ids]].astype('float64');l=local_bank()[d['bank_index'][ids]].astype('float64');c=d['cell'][ids].astype('int64');line=d['line'][ids].astype('int64')
    cf=np.eye(len(d['cell_names']))[c]*np.array(r['supported_cells'])[c,None];lf=np.eye(len(d['cell_line_names']))[line]*np.array(r['supported_lines'])[line,None]
    med,support,w,center,q,meta=transformations(d,train,True)
    assert r['activity_and_precision_transformations']==meta
    batch={k:d[k][ids].astype('float64') for k in ['effect','signed_effect','activity_reference','activity_alternate']};batch['activity_center']=center[ids].astype('float64');weights=w[ids].astype('float64');quality=q[ids].astype('float64')
    model=LineAwareEffectNetwork(r['cell_medians'],r['line_medians'],r['line_parent_cell'],'frozen').double().eval();model.load_state_dict(load_file(r['checkpoint']),strict=True)
    s={k:v.detach().numpy().copy() for k,v in model.state_dict().items()};expected=numpy_forward(s,g,l,c,cf,lf)
    tensors=[torch.from_numpy(v) for v in [g[:,0],g[:,1],l[:,0],l[:,1],c,line,cf,lf]];got=model(*tensors);errors={k:float(np.max(abs(got[k].detach().numpy()-expected[k]))) for k in expected};assert max(errors.values())<1e-7,errors
    tb={k:torch.from_numpy(v) for k,v in batch.items()};tw=torch.from_numpy(weights);tq=torch.from_numpy(quality);loss=activity_loss(got,tb,tw,tq,8,1,1);independent=numpy_loss(expected,batch,weights,quality,8,1,1);assert abs(float(loss.detach())-independent)<1e-8
    loss.backward();checks=[]
    for name in ['line_activity.weight','mag_line_offset','mag_hidden.weight','mag_out.weight','mag_cell_bias']:
        p=dict(model.named_parameters())[name];assert p.grad is not None and torch.isfinite(p.grad).all();index=int(p.grad.abs().argmax());analytic=float(p.grad.reshape(-1)[index]);assert abs(analytic)>1e-10
        original=float(p.detach().reshape(-1)[index]);eps=1e-6
        with torch.no_grad():
            p.reshape(-1)[index]=original+eps;plus=float(activity_loss(model(*tensors),tb,tw,tq,8,1,1))
            p.reshape(-1)[index]=original-eps;minus=float(activity_loss(model(*tensors),tb,tw,tq,8,1,1));p.reshape(-1)[index]=original
        finite=(plus-minus)/(2*eps);relative=abs(finite-analytic)/max(abs(finite),abs(analytic),1e-12);assert relative<1e-4,(name,relative);checks.append({'parameter':name,'relative_gradient_error':relative})
    # Masking an unknown context must remove any dependence on its raw index.
    unknown=list(tensors);unknown[-1]=torch.zeros_like(unknown[-1]);first=model(*unknown)['magnitude'].detach();unknown[5]=torch.zeros_like(unknown[5]);second=model(*unknown)['magnitude'].detach();assert torch.equal(first,second)
    assert sha(r['checkpoint'])==r['checkpoint_sha256'];assert_previous_sources()
    write_json(OUT/'first_fitted_network_independent_audit.json',{'passed':True,'checkpoint_sha256':r['checkpoint_sha256'],'NumPy_forward_errors':errors,'NumPy_loss_error':abs(float(loss.detach())-independent),'fitted_gradient_checks':checks,'unknown_line_index_has_no_effect_when_context_is_masked':True,'actual_measured_examples_include':['Known signed Kosicki contrasts','New Retallick measured allele activities','New McAfee nearby and cross-cell effects'],'new_fitting_or_saved_weight_change':False,'scope':'Independent numerical and learned-gradient correctness only; no successful external prediction or human genetic validation is established.'})
    print('LINE-AWARE FITTED MODEL INDEPENDENTLY VERIFIED',flush=True)


if __name__=='__main__':main()

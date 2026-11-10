"""Independent loss/batch checks and a sustained real-data memory preflight.

Engineering updates are discarded. Protected outside-study outcomes are unread.
"""
import os
os.environ.setdefault('PYTORCH_MPS_HIGH_WATERMARK_RATIO','0.4')
os.environ.setdefault('PYTORCH_MPS_LOW_WATERMARK_RATIO','0.3')
import ctypes
import shutil
import struct
import time
from collections import Counter
from bounded_common import *
from native_data import load_data,PackedDNA,normalization
from bounded_network import NativeAdaptationModel,torch
from native_network import NativeAdaptationModel as OriginalModel
from stable_loss import activity_loss
from activity_effect_network import activity_loss as original_loss
from stable_batches import epoch_batches,pair_buckets,assemble
from training_utilities import transformations


def memory(device='cpu'):
    # Darwin rusage_info_v2: UUID then 64-bit counters; physical footprint is
    # counter 7. The kernel's footprint includes compressed resident pages.
    lib=ctypes.CDLL('/usr/lib/libproc.dylib');buf=ctypes.create_string_buffer(1024)
    code=lib.proc_pid_rusage(os.getpid(),2,ctypes.byref(buf));assert code==0
    values=struct.unpack_from('10Q',buf.raw,16)
    result={'physical_footprint_bytes':values[7],
            'free_disk_bytes':shutil.disk_usage(ROOT).free}
    if device=='mps':
        torch.mps.synchronize()
        result.update(mps_current_bytes=torch.mps.current_allocated_memory(),
                      mps_driver_bytes=torch.mps.driver_allocated_memory())
    return result


def independent_loss(pred,batch,w,q,nmain,nlocal,ncross):
    """NumPy reference, deliberately retaining original subset reductions."""
    p={k:v.detach().numpy() for k,v in pred.items()}
    b={k:v.detach().numpy() for k,v in batch.items()};w=w.numpy();q=q.numpy()
    result=np.mean(w[:nmain]*np.abs(p['magnitude'][:nmain]-b['effect'][:nmain]))
    signed=np.flatnonzero(np.isfinite(b['signed_effect'][:nmain]))
    if len(signed):result+=.25*np.mean(w[signed]*q[signed]*np.abs(p['signed_effect'][signed]-b['signed_effect'][signed]))
    active=np.flatnonzero(np.isfinite(b['activity_reference'][:nmain])&np.isfinite(b['activity_alternate'][:nmain]))
    def huber(v):return np.where(np.abs(v)<1,.5*v*v,np.abs(v)-.5)
    if len(active):
        a=b['activity_reference'][active]-b['activity_center'][active]
        z=b['activity_alternate'][active]-b['activity_center'][active]
        result+=.05*np.mean(w[active]*(huber(p['activity_reference'][active]-a)+huber(p['activity_alternate'][active]-z))*.5)
        result+=.05*np.mean(w[active]*q[active]*huber(p['signed_effect'][active]-(z-a)))
    for n,offset in [(nlocal,nmain),(ncross,nmain+2*nlocal)]:
        if n:
            pp=p['magnitude'][offset:offset+2*n].reshape(n,2);y=b['effect'][offset:offset+2*n].reshape(n,2)
            qq=np.sqrt(q[offset:offset+2*n].reshape(n,2).prod(1))
            result+=.2*np.mean(qq*np.abs((pp[:,0]-pp[:,1])-(y[:,0]-y[:,1])))
    return float(result)


def main():
    start=time.monotonic();assert_versions_stable();plan=get_stable_plan()
    out=DEST/'stable_runtime_preflight.json';assert not out.exists(),'Inspect a completed preflight before repeating it.'
    d,pairs,bank,x,extra=load_data();first,neural=add_all_documented_signed_targets(d)
    packed=PackedDNA(d,bank);train=d['split']==0;tr=np.flatnonzero(train)
    med,support,w,centers,quality,context=transformations(d,train,True)
    lengths=packed.length[d['bank_index']];opt=plan['optimization'];pools=pair_buckets(pairs,lengths,train)
    nl=len(d['cell_line_names']);parent=np.zeros(nl,dtype='int64');lm=np.zeros(nl,dtype='float32');ls=np.zeros(nl,dtype=bool)
    for line in range(nl):
        cells=np.unique(d['cell'][d['line']==line]);assert len(cells)<=1
        if len(cells):parent[line]=cells[0]
        good=train&(d['line']==line);ls[line]=good.any();lm[line]=np.median(d['effect'][good]) if good.any() else med[parent[line]]
    norm=normalization(bank,x,extra,np.unique(d['bank_index'][tr]));rng=np.random.default_rng(42)
    batches=epoch_batches(tr,lengths,rng,opt['main_batch']);replay=np.concatenate([ids for length,ids in batches])
    assert np.array_equal(np.sort(replay),tr)
    rng2=np.random.default_rng(42);other=epoch_batches(tr,lengths,rng2,opt['main_batch'])
    assert all(l==ll and np.array_equal(a,b) for (l,a),(ll,b) in zip(batches,other))
    concrete=[];shape_examples={};schedule1=hashlib.sha256();schedule2=hashlib.sha256()
    for (length,main),(length2,main2) in zip(batches,other):
        ids,counts=assemble(main,length,pools,rng,opt['local_pairs'],opt['cross_pairs'])
        ids2,counts2=assemble(main2,length2,pools,rng2,opt['local_pairs'],opt['cross_pairs'])
        assert counts==counts2 and np.array_equal(ids,ids2)
        assert train[ids].all() and np.all(lengths[ids]==length)
        schedule1.update(ids.tobytes());schedule2.update(ids2.tobytes())
        entry=(length,main,ids,counts);concrete.append(entry)
        shape_examples.setdefault((length,len(main),*counts),entry)
    assert schedule1.hexdigest()==schedule2.hexdigest()
    # Each validation row has one unique output slot after length sorting.
    te=np.flatnonzero(d['split']==1);inverse=np.full(len(lengths),-1,dtype='int32');inverse[te]=np.arange(len(te))
    visited=np.concatenate([te[lengths[te]==length] for length in np.unique(lengths[te])])
    assert np.array_equal(np.sort(inverse[visited]),np.arange(len(te)))
    # Actual-data losses plus deliberately empty supervision masks. Synthetic
    # predictions here test algebra, never serve as biological training labels.
    examples=np.r_[tr[neural[tr]][:16],tr[first[tr]][:16],tr[:32]]
    nmain=48;nlocal=4;ncross=4;assert len(examples)==64
    target={k:torch.from_numpy(d[k][examples]).double() for k in ['effect','signed_effect','activity_reference','activity_alternate']}
    target['activity_center']=torch.from_numpy(centers[examples]).double()
    tw=torch.from_numpy(w[examples]).double();tq=torch.from_numpy(quality[examples]).double()
    torch.manual_seed(725);base={k:torch.randn(64,dtype=torch.float64) for k in ['magnitude','signed_effect','activity_reference','activity_alternate']}
    loss_checks=[]
    for mode in ['actual','no_signed','no_activity','neither','all_known']:
        b={k:v.clone() for k,v in target.items()}
        if mode in ['no_signed','neither']:b['signed_effect'][:]=float('nan')
        if mode in ['no_activity','neither']:
            b['activity_reference'][:]=float('nan');b['activity_alternate'][:]=float('nan')
        if mode=='all_known':
            for k in ['signed_effect','activity_reference','activity_alternate']:b[k]=torch.nan_to_num(b[k],nan=.3)
        p={k:v.clone().requires_grad_() for k,v in base.items()};p2={k:v.clone().requires_grad_() for k,v in base.items()}
        old=original_loss(p,b,tw,tq,nmain,nlocal,ncross);new=activity_loss(p2,b,tw,tq,nmain,nlocal,ncross)
        ref=independent_loss(p,b,tw,tq,nmain,nlocal,ncross)
        old.backward();new.backward();errors=[]
        for k in p:
            a=p[k].grad if p[k].grad is not None else torch.zeros_like(p[k]);z=p2[k].grad if p2[k].grad is not None else torch.zeros_like(p2[k])
            errors.append(float((a-z).abs().max()));assert torch.isfinite(z).all()
        error=max(abs(float(old.detach())-ref),abs(float(new.detach())-ref),*errors);assert error<1e-10,(mode,error)
        loss_checks.append({'mask_case':mode,'maximum_loss_or_gradient_error':error})
    torch.manual_seed(42);old=OriginalModel(med,lm,parent,norm,'adapt_last_stage').double().eval()
    torch.manual_seed(42);new=NativeAdaptationModel(med,lm,parent,norm,'adapt_last_stage').double().eval()
    assert all(torch.equal(v,new.state_dict()[k]) for k,v in old.state_dict().items())
    chosen=np.r_[tr[neural[tr]][:4],tr[neural[tr]][:4],*[tr[lengths[tr]==l][:2] for l in np.unique(lengths[tr])]]
    with torch.no_grad():
        a=old(chosen,d,packed,bank,x,extra,support,ls);b=new(chosen,d,packed,bank,x,extra,support,ls)
    forward_error=max(float((a[k]-b[k]).abs().max()) for k in a);assert forward_error<1e-10
    # Independent original full-batch versus bounded-reader gradients on
    # true neuron measurements. Both use double precision, fixed BN and no
    # head dropout. Extra unscored copies must have exactly zero influence.
    gradient_ids=np.r_[tr[neural[tr]][:3],tr[neural[tr]][:2]]
    gb={k:torch.from_numpy(d[k][gradient_ids]).double() for k in ['effect','signed_effect','activity_reference','activity_alternate']}
    gb['activity_center']=torch.from_numpy(centers[gradient_ids]).double()
    gw=torch.from_numpy(w[gradient_ids]).double();gq=torch.from_numpy(quality[gradient_ids]).double()
    old.train();new.train();old.zero_grad();new.zero_grad();torch.manual_seed(777)
    original_loss(old(gradient_ids,d,packed,bank,x,extra,support,ls),gb,gw,gq,len(gradient_ids),0,0).backward()
    torch.manual_seed(777)
    activity_loss(new(gradient_ids,d,packed,bank,x,extra,support,ls),gb,gw,gq,len(gradient_ids),0,0).backward()
    gradient_checks=[]
    for name,parameter in old.named_parameters():
        if not parameter.requires_grad:continue
        other=dict(new.named_parameters())[name]
        ga=parameter.grad if parameter.grad is not None else torch.zeros_like(parameter)
        gz=other.grad if other.grad is not None else torch.zeros_like(other)
        error=float((ga-gz).abs().max());scale=max(float(ga.abs().max()),float(gz.abs().max()),1e-8)
        assert error<1e-8+1e-7*scale,(name,error,scale)
        gradient_checks.append({'parameter':name,'maximum_absolute_error':error,'gradient_scale':scale})
    write_json(DEST/'bounded_calculation_equivalence.json',{'passed':True,'loss_checks':loss_checks,
        'duplicate_batch_forward_max_error':forward_error,'original_vs_bounded_gradient_checks':gradient_checks,
        'unscored_real_construct_copies_do_not_change_supervision':True,'new_test_outcomes_read':False})
    print('BOUNDED NATIVE CALCULATION EQUIVALENCE PASS',flush=True)
    del old,new,a,b,p,p2
    import gc;gc.collect()
    device=plan['resource_revision']['device']
    # Worst-case adapted arm, all real length/tail/pair shapes covered first;
    # remaining steps replay fixed scheduled batches, not selected outcomes.
    torch.manual_seed(42);model=NativeAdaptationModel(med,lm,parent,norm,'adapt_last_stage').to(device).train()
    model.anchored={k:v.to(device) for k,v in model.anchored.items()}
    initial={k:v.detach().cpu().clone() for k,v in model.reader.state_dict().items()}
    optimizer=torch.optim.AdamW([{'params':list(model.head.parameters()),'lr':opt['head_lr']},
        {'params':[p for p in model.reader.parameters() if p.requires_grad],'lr':opt['reader_lr']}],weight_decay=opt['weight_decay'])
    order=list(shape_examples.values())+concrete;observations=[]
    for step in range(1000):
        length,main,ids,counts=order[step%len(order)]
        optimizer.zero_grad(set_to_none=True);pred=model(ids,d,packed,bank,x,extra,support,ls)
        batch={k:torch.from_numpy(d[k][ids]).to(device) for k in ['effect','signed_effect','activity_reference','activity_alternate']}
        batch['activity_center']=torch.from_numpy(centers[ids]).to(device)
        loss=activity_loss(pred,batch,torch.from_numpy(w[ids]).to(device),torch.from_numpy(quality[ids]).to(device),len(main),*counts)+model.anchor_loss()
        assert torch.isfinite(loss);loss.backward()
        assert all(p.grad is None or torch.isfinite(p.grad).all() for p in model.parameters())
        torch.nn.utils.clip_grad_norm_([p for p in model.parameters() if p.requires_grad],opt['gradient_clip']);optimizer.step()
        if (step+1)%100==0 and device=='mps':torch.mps.synchronize();torch.mps.empty_cache()
        if step==0 or (step+1)%10==0:
            row={'step':step+1,'elapsed_seconds':time.monotonic()-start,**memory(device)};observations.append(row)
            write_json(DEST/'stable_runtime_preflight_live.json',{'PID':os.getpid(),'observations':observations,'passed':False,'status':'running'})
            print('STABLE PREFLIGHT',json.dumps(row),flush=True)
            assert row['physical_footprint_bytes']<3*1024**3,'Resource ceiling exceeded; preserve observations.'
            assert row['free_disk_bytes']>256*1024**2,'Insufficient disk reserve; preserve observations.'
    # Every validation tensor shape is exercised and original slots restored.
    model.eval();restored=np.full(len(te),np.nan,dtype='float32')
    with torch.no_grad():
        for length in np.unique(lengths[te]):
            selected=te[lengths[te]==length]
            for i in range(0,len(selected),256):
                ix=selected[i:i+256];p=model(ix,d,packed,bank,x,extra,support,ls)
                restored[inverse[ix]]=p['magnitude'].cpu().numpy()
    assert np.isfinite(restored).all()
    after=model(chosen,d,packed,bank,x,extra,support,ls);reverse=model(chosen,d,packed,bank,x,extra,support,ls,True)
    swap_error=max(float((after['magnitude']-reverse['magnitude']).abs().max().detach().cpu()),
                   float((after['signed_effect']+reverse['signed_effect']).abs().max().detach().cpu()))
    assert swap_error<1e-5
    changes={k:int(torch.count_nonzero(v.detach().cpu()!=initial[k])) for k,v in model.reader.state_dict().items()}
    assert sum(v for k,v in changes.items() if k.startswith('main.3.'))>0
    assert sum(v for k,v in changes.items() if not k.startswith('main.3.'))==0
    assert sum(v for k,v in changes.items() if 'running_' in k or 'num_batches' in k)==0
    # Check plateau after shape-cache warmup, not only a short three-step run.
    late_growth=max(r['physical_footprint_bytes'] for r in observations if r['step']>=750)-max(r['physical_footprint_bytes'] for r in observations if 250<=r['step']<=500)
    final=memory(device);criteria={'thousand_real_backprop_updates':True,'late_growth_under_256_MiB':late_growth<256*1024**2,
        'peak_footprint_under_3_GiB':max(r['physical_footprint_bytes'] for r in observations)<3*1024**3,
        'final_footprint_under_3_GiB':final['physical_footprint_bytes']<3*1024**3,'disk_reserve_over_256_MiB':final['free_disk_bytes']>256*1024**2}
    receipt={'passed':all(criteria.values()),'criteria':criteria,'device':device,'PID':os.getpid(),
       'loss_and_prediction_gradient_equivalence':loss_checks,'duplicate_batch_forward_max_error':forward_error,
       'every_main_training_row_exactly_once':True,'training_rows':len(tr),'batch_schedule_identical_between_arms':True,
       'schedule_sha256':schedule1.hexdigest(),'bounded_training_shapes':[list(k) for k in shape_examples],
       'all_validation_rows_restored':True,'observations':observations,'late_footprint_growth_bytes':late_growth,'final_memory':final,
       'allele_swap_max_error':swap_error,'reader_changed_scalar_counts':changes,
       'engineering_updates_discarded':True,'new_outside_study_outcomes_read':False,
       'scientific_prediction_success_established':False,'elapsed_seconds':time.monotonic()-start,
       'code_sha256':{p.name:sha(p) for p in Path(__file__).parent.glob('*.py')},'original_vs_bounded_gradient_checks':gradient_checks}
    assert_versions_stable();write_json(out,receipt);print('STABLE RUNTIME PREFLIGHT FINISHED',receipt['passed'],flush=True)
    assert receipt['passed']


if __name__=='__main__':main()

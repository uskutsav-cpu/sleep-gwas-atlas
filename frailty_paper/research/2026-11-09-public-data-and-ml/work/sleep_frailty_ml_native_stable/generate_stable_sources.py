"""Preserve the original trial and derive only the execution-resource revision."""
from stable_common import *


def replace_once(code,old,new):
    assert code.count(old)==1,(old,code.count(old));return code.replace(old,new)


def main():
    assert not (DEST/'cpu_resource_revision_receipt.json').exists(), 'CPU revision already exists; this historical MPS generator must not overwrite it.'
    assert_versions_stable();get_stable_plan();source=ROOT/'work/sleep_frailty_ml_native';dest=Path(__file__).parent;changes={}
    code=(source/'native_network.py').read_text().replace('from native_common import *','from stable_common import *')
    code=replace_once(code,'unique,inverse=np.unique(b,return_inverse=True)','unique=b;inverse=np.arange(len(b),dtype="int64")')
    (dest/'stable_network.py').write_text(code);changes['stable_network.py']='Identical network and features; evaluate actual repeated constructs separately so true-length batches retain fixed shapes.'
    code=(source/'train_native_development.py').read_text()
    code=code.replace('from native_common import *','from stable_common import *').replace('from native_network import','from stable_network import').replace('from activity_effect_network import activity_loss','from stable_loss import activity_loss\nfrom stable_batches import epoch_batches,pair_buckets,assemble')
    code=code.replace('assert_versions_native()','assert_versions_stable()').replace('get_plan()','get_stable_plan()')
    code=code.replace("DEST/'native_runtime_and_gradient_audit.json'","DEST/'stable_runtime_preflight.json'").replace("DEST/'finite_native_adaptation_plan.json'","DEST/'finite_stable_native_plan.json'")
    code=code.replace("Path(__file__).with_name('native_network.py')","Path(__file__).with_name('stable_network.py')").replace("Path(__file__).with_name('native_data.py')","ROOT/'work/sleep_frailty_ml_native/native_data.py'")
    old="optimizer=torch.optim.AdamW(groups,weight_decay=opt['weight_decay']);tp={k:p[np.all(train[p],axis=1)] for k,p in pairs.items()}"
    new="optimizer=torch.optim.AdamW(groups,weight_decay=opt['weight_decay']);lengths=packed.length[d['bank_index']];tp=pair_buckets(pairs,lengths,train)"
    code=replace_once(code,old,new)
    old="""            for i in range(0,len(te),256):
                p=model(te[i:i+256],d,packed,bank,x,extra,support,ls)
                for k in result:result[k][i:i+256]=p[k].cpu().numpy()"""
    new="""            inverse=np.full(len(d['effect']),-1,dtype='int32');inverse[te]=np.arange(len(te))
            for length in np.unique(lengths[te]):
                selected=te[lengths[te]==length]
                for i in range(0,len(selected),256):
                    ix=selected[i:i+256];p=model(ix,d,packed,bank,x,extra,support,ls)
                    for k in result:result[k][inverse[ix]]=p[k].cpu().numpy()"""
    code=replace_once(code,old,new)
    old="""        model.train();order=rng.permutation(tr);objectives=[];anchors=[]
        for i in range(0,len(order),opt['main_batch']):
            main_ids=order[i:i+opt['main_batch']];chunks=[main_ids];counts=[]
            for key,n in [('local',opt['local_pairs']),('cross',opt['cross_pairs'])]:
                p=tp[key];count=n if len(p) else 0;counts.append(count)
                if count:chunks.append(p[rng.integers(len(p),size=count)].reshape(-1))
            ids=np.concatenate(chunks);schedule.update(ids.tobytes());optimizer.zero_grad();prediction=model(ids,d,packed,bank,x,extra,support,ls)"""
    new="""        model.train();batches=epoch_batches(tr,lengths,rng,opt['main_batch']);objectives=[];anchors=[];completed_rows=0
        replay=np.concatenate([rows for length,rows in batches]);assert np.array_equal(np.sort(replay),tr)
        for length,main_ids in batches:
            ids,counts=assemble(main_ids,length,tp,rng,opt['local_pairs'],opt['cross_pairs']);assert np.all(lengths[ids]==length)
            completed_rows+=len(main_ids);schedule.update(ids.tobytes());optimizer.zero_grad();prediction=model(ids,d,packed,bank,x,extra,support,ls)"""
    code=replace_once(code,old,new)
    code=code.replace("'epoch_rows_completed':min(i+len(main_ids),len(tr))","'epoch_rows_completed':completed_rows")
    old="objectives.append(float(objective.detach().cpu()));anchors.append(float(anchor.detach().cpu()))"
    new=old+"\n            if device=='mps' and updates%100==0:torch.mps.empty_cache()"
    code=replace_once(code,old,new)
    code=replace_once(code,'        p=predict();metrics=',"        if device=='mps':torch.mps.synchronize();torch.mps.empty_cache()\n        p=predict();metrics=")
    code="import os\nos.environ.setdefault('PYTORCH_MPS_HIGH_WATERMARK_RATIO','0.4')\nos.environ.setdefault('PYTORCH_MPS_LOW_WATERMARK_RATIO','0.3')\n"+code
    (dest/'train_stable_development.py').write_text(code);changes['train_stable_development.py']='True-length grouped shuffled batches and same-length pairs; complete-row coverage verified; fixed-shape masked objective; true-length prediction restored to original order; bounded MPS cache policy.'
    code=(source/'run_native_comparison.py').read_text().replace('from native_common import *','from stable_common import *').replace('assert_versions_native()','assert_versions_stable()').replace('get_plan()','get_stable_plan()')
    code=code.replace("'train_native_development.py'","'train_stable_development.py'").replace("DEST/'native_runtime_and_gradient_audit.json'","DEST/'stable_runtime_preflight.json'")
    code=code.replace("'VECLIB_MAXIMUM_THREADS':'2'","'VECLIB_MAXIMUM_THREADS':'2','PYTORCH_MPS_HIGH_WATERMARK_RATIO':'0.4','PYTORCH_MPS_LOW_WATERMARK_RATIO':'0.3'")
    (dest/'run_stable_comparison.py').write_text(code);changes['run_stable_comparison.py']='Same neural gate, phase roles, three epochs, seeds and protected new test, in a new version with fresh original initialization.'
    write_json(DEST/'stable_source_derivation_receipt.json',{'original_code_sha256':{p.name:sha(p) for p in source.glob('*.py')},
       'derived_code_sha256':{p.name:sha(p) for p in dest.glob('*.py')},'changes':changes,'old_partial_run_and_first_checkpoint_preserved':True})
    assert_versions_stable();print('BOUNDED-SHAPE RESOURCE REVISION SOURCES WRITTEN',flush=True)


if __name__=='__main__':main()

import os
os.environ.setdefault('PYTORCH_MPS_HIGH_WATERMARK_RATIO','0.4')
os.environ.setdefault('PYTORCH_MPS_LOW_WATERMARK_RATIO','0.3')
"""Separately gated, matched permissible-score input training; not yet executed."""
import argparse,time,shutil
from resumption import save_recovery,load_recovery
from resumable_common import *
from native_data import PackedDNA,normalization
from readonly_native_inputs import load_readonly_inputs
from permitted_AVI_features import fit_scaling
from permitted_AVI_experiment import get_AVI_plan,PLAN_PATH
from permitted_AVI_network import PermittedAVIModel,torch
from training_utilities import transformations,detailed_metrics
from stable_loss import activity_loss
from stable_batches import epoch_batches,pair_buckets,assemble
from safetensors.torch import save_file


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--recipe',required=True);parser.add_argument('--held-study',required=True);parser.add_argument('--resume',action='store_true');parser.add_argument('--arm',choices=['scoreless','scores'],required=True)
    parser.add_argument('--audit-tag');parser.add_argument('--audit-stop-after-updates',type=int);args=parser.parse_args()
    if args.audit_tag:assert args.audit_tag.replace('_','').isalnum() and len(args.audit_tag)<80 and args.audit_stop_after_updates is not None
    if args.audit_stop_after_updates is not None:assert args.audit_tag and args.audit_stop_after_updates>0
    assert_versions_stable();plan=get_AVI_plan(for_execution_audit=bool(args.audit_tag));assert args.recipe in plan['recipes'] and args.held_study in plan['phase_one_roles']+plan['phase_two_roles_if_gate_passes']
    proof=json.loads((BOUNDED/'stable_runtime_preflight.json').read_text());assert proof['passed']
    if not args.audit_tag:assert json.loads((DEST/'recovery_state_audit.json').read_text())['passed']
    for filename in ['train_stable_development.py','bounded_network.py','bounded_common.py']:
        assert sha(ROOT/'work/sleep_frailty_ml_bounded'/filename)==proof['code_sha256'][filename]
    assert shutil.disk_usage(ROOT).free>=plan['launch_system_reserve_bytes'],'System storage below the separate AVI launch reserve; scientific fit not started.'
    d,pairs,bank,x,extra,input_receipt=load_readonly_inputs();added_siraj,added_neural=add_all_documented_signed_targets(d);packed=PackedDNA(d,bank);train=d['split']==0
    if args.held_study=='none':test=d['split']==1
    else:
        sid=list(d['study_names']).index(args.held_study);groups=np.unique(d['group'][d['study']==sid]);train&=(d['study']!=sid)&~np.isin(d['group'],groups);test=(d['study']==sid)&(d['split']!=2)
    tr,te=np.flatnonzero(train),np.flatnonzero(test);assert not set(d['group'][tr])&set(d['group'][te])
    med,support,w,centers,quality,context=transformations(d,train,True);norm=normalization(bank,x,extra,np.unique(d['bank_index'][tr]))
    nl=len(d['cell_line_names']);parent=np.zeros(nl,dtype='int64');lm=np.zeros(nl,dtype='float32');ls=np.zeros(nl,dtype=bool)
    for line in range(nl):
        cells=np.unique(d['cell'][d['line']==line]);assert len(cells)<=1
        if len(cells):parent[line]=cells[0]
        ix=train&(d['line']==line);ls[line]=ix.any();lm[line]=np.median(d['effect'][ix]) if ix.any() else med[parent[line]]
    avi_path=DEST/'AVI_preparation/corrected_permitted_AVI_physical_bank.npy';assert sha(avi_path)==plan['permitted_score_bank_sha256'];avi_bank=np.load(avi_path,mmap_mode='r')
    avi_scaling=fit_scaling(avi_bank,d['bank_index'][tr]);scaling_proof=json.loads((DEST/'AVI_preparation/permitted_AVI_train_only_scaling_audit.json').read_text());assert avi_scaling==scaling_proof['roles'][args.held_study]['scaling']
    seed=plan['seed'];torch.manual_seed(seed);device=plan['resource_revision']['device']
    model=PermittedAVIModel(med,lm,parent,norm,args.recipe,avi_scaling,args.arm=='scores',seed).to(device)
    model.anchored={k:v.to(device) for k,v in model.anchored.items()}
    initial_reader={k:v.detach().cpu().clone() for k,v in model.reader.state_dict().items()};h=hashlib.sha256()
    for key,value in sorted(model.state_dict().items()):h.update(key.encode());h.update(value.detach().cpu().numpy().tobytes())
    directory=(DEST/'AVI_recovery_audit_runs'/args.audit_tag/'development') if args.audit_tag else DEST/'AVI_development'
    directory.mkdir(exist_ok=True,parents=True);name=f'{args.recipe}_{args.held_study}_{args.arm}_{seed}';record_path=directory/(name+'_record.json');assert not record_path.exists()
    assert (directory/(name+'_started.json')).exists()==args.resume,'A fresh fit and a resumed fit require distinct explicit state.'
    meta={'recipe':args.recipe,'arm':args.arm,'AVI_scaling':avi_scaling,'AVI_bank_sha256':sha(avi_path),'readonly_input_receipt':input_receipt,'held_study':args.held_study,'seed':seed,'epochs':plan['epochs'],'training_rows':len(tr),'prediction_rows':len(te),
       'training_index_sha256':hashlib.sha256(tr.tobytes()).hexdigest(),'prediction_index_sha256':hashlib.sha256(te.tobytes()).hexdigest(),
       'initial_parameter_sha256':h.hexdigest(),'normalization':norm,'cell_medians':med.tolist(),'line_medians':lm.tolist(),'line_parent':parent.tolist(),
       'supported_cells':support.tolist(),'supported_lines':ls.tolist(),'context_transformations':context,
       'new_original_Siraj_signed_training_rows':int(added_siraj[tr].sum()),'new_original_neuronal_signed_training_rows':int(added_neural[tr].sum()),
       'total_known_signed_training_rows':int(np.isfinite(d['signed_effect'][tr]).sum()),
       'trainable_head_parameters':sum(p.numel() for p in model.head.parameters()),
       'trainable_reader_parameters':sum(p.numel() for p in model.reader.parameters() if p.requires_grad),
       'plan_sha256':sha(PLAN_PATH),'pretrained_reader_sha256':plan['reader_sha256'],
       'code_sha256':{p.name:sha(p) for p in [Path(__file__),Path(__file__).with_name('permitted_AVI_network.py'),Path(__file__).with_name('permitted_AVI_head.py'),Path(__file__).with_name('permitted_AVI_features.py'),Path(__file__).with_name('readonly_native_inputs.py'),Path(__file__).with_name('permitted_AVI_experiment.py'),ROOT/'work/sleep_frailty_ml_bounded/bounded_network.py',ROOT/'work/sleep_frailty_ml_native/native_data.py']},
       'device':device,'torch_version':torch.__version__,'physical_source_receipt':packed.receipt,
       'source_standard_error_used_at_inference':False,'calibration_rows_used_for_fit':False,'reserved_Ulirsch_outcomes_used':False,'PID':os.getpid()}
    identity={k:meta[k] for k in ['recipe','arm','AVI_scaling','AVI_bank_sha256','readonly_input_receipt','held_study','seed','epochs','training_index_sha256','prediction_index_sha256','initial_parameter_sha256','plan_sha256','code_sha256']}
    identity['recovery_helper_sha256']=sha(Path(__file__).with_name('resumption.py'))
    if args.resume:
        original=json.loads((directory/(name+'_started.json')).read_text())
        for key in identity:
            if key!='recovery_helper_sha256':assert original[key]==meta[key],key
    else:write_json(directory/(name+'_started.json'),meta)
    opt=plan['optimization']
    groups=[{'params':list(model.head.parameters()),'lr':opt['head_lr']}]
    core=[p for p in model.reader.parameters() if p.requires_grad]
    if core:groups.append({'params':core,'lr':opt['reader_lr']})
    optimizer=torch.optim.AdamW(groups,weight_decay=opt['weight_decay']);lengths=packed.length[d['bank_index']];tp=pair_buckets(pairs,lengths,train)
    baseline=np.where(ls[d['line'][te]],lm[d['line'][te]],med[d['cell'][te]])
    model.eval()
    with torch.no_grad():assert abs(model(te[:32],d,packed,bank,x,extra,support,ls,avi_bank)['magnitude'].cpu().numpy()-baseline[:32]).max()<1e-6
    def predict():
        model.eval();result={k:np.empty(len(te),dtype='float32') for k in ['magnitude','signed_effect','activity_reference','activity_alternate']}
        with torch.no_grad():
            inverse=np.full(len(d['effect']),-1,dtype='int32');inverse[te]=np.arange(len(te))
            for length in np.unique(lengths[te]):
                selected=te[lengths[te]==length]
                for i in range(0,len(selected),256):
                    ix=selected[i:i+256];p=model(ix,d,packed,bank,x,extra,support,ls,avi_bank)
                    for k in result:result[k][inverse[ix]]=p[k].cpu().numpy()
        return result
    recovery=directory/(name+'_recovery.pt')
    rng=np.random.default_rng(seed);schedule=hashlib.sha256();updates=0;curve=[];start=time.monotonic();last_print=start
    restored=load_recovery(recovery,model,optimizer,device,identity) if args.resume else None
    saved_updates=int(restored['updates']) if restored else 0;updates=saved_updates
    saved_epoch=int(restored['epoch']) if restored else 0
    elapsed_before=float(restored['elapsed_seconds']) if restored else 0.
    curve=restored['curve'] if restored else []
    if restored:
        # Preserve any writes made after the last atomically committed state.
        beyond=[]
        for path in directory.glob(name+'_epoch*'):
            number=int(path.name.split('_epoch',1)[1].split('_',1)[0].split('.',1)[0])
            if number>len(curve):beyond.append(path)
        curvepath=directory/(name+'_curve.json')
        if curvepath.exists():beyond.append(curvepath)
        if beyond:
            archive=directory/('recovery_archive_'+str(time.time_ns()));archive.mkdir()
            manifest={p.name:sha(p) for p in beyond}
            for path in beyond:path.rename(archive/path.name)
            write_json(archive/'preserved_after_last_committed_state.json',manifest)
        history=directory/(name+'_resume_history.json')
        events=json.loads(history.read_text()) if history.exists() else []
        events.append({'PID':os.getpid(),'recovery_sha256':sha(recovery),'updates':saved_updates,'epoch':saved_epoch,
            'model_optimizer_and_CPU_MPS_random_state_restored_bit_for_bit':restored['restore_verified_exactly']})
        write_json(history,events)
        if curve:write_json(curvepath,curve)
    ordinal=0;chain_verified=saved_updates==0
    if restored and saved_updates==0:
        assert restored['schedule_sha256']==schedule.hexdigest()
        assert json.loads(restored['NumPy_rng_json'])==rng.bit_generator.state
    epoch=1;objectives=[];anchors=[]
    def checkpoint_now(phase):
        save_recovery(recovery,model,optimizer,device,identity,{'updates':updates,'epoch':epoch,'phase':phase,
           'schedule_sha256':schedule.hexdigest(),'NumPy_rng_json':json.dumps(rng.bit_generator.state),
           'objectives':objectives,'anchors':anchors,'curve':curve,'elapsed_seconds':elapsed_before+time.monotonic()-start})
        assert shutil.disk_usage(ROOT).free>300*1024**2,'Stop with committed recovery state before exhausting disk reserve.'
    if not restored:
        torch.manual_seed(seed+1)
        if device=='mps':torch.mps.manual_seed(seed+1)
        checkpoint_now('training')
    for epoch in range(1,plan['epochs']+1):
        if not restored or epoch>saved_epoch:
            torch.manual_seed(seed+epoch)
            if device=='mps':torch.mps.manual_seed(seed+epoch)
        model.train();batches=epoch_batches(tr,lengths,rng,opt['main_batch']);completed_rows=0
        objectives=list(restored['objectives']) if restored and epoch==saved_epoch and restored['phase']=='training' else []
        anchors=list(restored['anchors']) if restored and epoch==saved_epoch and restored['phase']=='training' else []
        replay=np.concatenate([rows for length,rows in batches]);assert np.array_equal(np.sort(replay),tr)
        for length,main_ids in batches:
            ids,counts=assemble(main_ids,length,tp,rng,opt['local_pairs'],opt['cross_pairs']);assert np.all(lengths[ids]==length)
            completed_rows+=len(main_ids);schedule.update(ids.tobytes());ordinal+=1
            if ordinal<=saved_updates:
                if ordinal==saved_updates:
                    assert epoch==saved_epoch
                    assert schedule.hexdigest()==restored['schedule_sha256']
                    assert rng.bit_generator.state==json.loads(restored['NumPy_rng_json'])
                    chain_verified=True
                continue
            assert chain_verified and updates==ordinal-1
            optimizer.zero_grad();prediction=model(ids,d,packed,bank,x,extra,support,ls,avi_bank)
            batch={k:torch.from_numpy(d[k][ids]).to(device) for k in ['effect','signed_effect','activity_reference','activity_alternate']};batch['activity_center']=torch.from_numpy(centers[ids]).to(device)
            objective=activity_loss(prediction,batch,torch.from_numpy(w[ids]).to(device),torch.from_numpy(quality[ids]).to(device),len(main_ids),*counts)
            anchor=model.anchor_loss();loss=objective+anchor;assert torch.isfinite(loss);loss.backward()
            assert all(p.grad is None or torch.isfinite(p.grad).all() for p in model.parameters())
            torch.nn.utils.clip_grad_norm_([p for p in model.parameters() if p.requires_grad],opt['gradient_clip']);optimizer.step();updates+=1
            objectives.append(float(objective.detach().cpu()));anchors.append(float(anchor.detach().cpu()))
            if args.audit_stop_after_updates and updates>=args.audit_stop_after_updates:
                checkpoint_now('training')
                write_json(directory/(name+'_audit_stop.json'),{'audit_only':True,'updates':updates,'schedule_sha256':schedule.hexdigest(),'used_to_initialize_scientific_fit':False})
                print('RECOVERY AUDIT SEGMENT FINISHED',args.audit_tag,updates,flush=True)
                return
            if device=='mps' and updates%100==0:torch.mps.empty_cache()
            if updates%plan['recovery_execution']['checkpoint_every_updates']==0:checkpoint_now('training')
            if time.monotonic()-last_print>30:
                assert shutil.disk_usage(ROOT).free>300*1024**2,'Disk reserve exhausted during fit; preserve partial results and inspect before restarting.'
                progress={'fit':name,'epoch':epoch,'epoch_rows_completed':completed_rows,
                   'epoch_rows_total':len(tr),'updates':updates,'elapsed_seconds':elapsed_before+time.monotonic()-start,'PID':os.getpid()}
                write_json(directory/(name+'_live_progress.json'),progress)
                print('NATIVE FIT PROGRESS',json.dumps(progress),flush=True);last_print=time.monotonic()
        if restored and (epoch<saved_epoch or (epoch==saved_epoch and restored['phase']=='epoch_complete')):continue
        if device=='mps':torch.mps.synchronize();torch.mps.empty_cache()
        p=predict();assert shutil.disk_usage(ROOT).free>300*1024**2,'Insufficient disk reserve for epoch artifacts.'
        metrics=detailed_metrics(d,te,p,med,support,pairs);skills=[]
        for c in np.unique(d['cell'][te]):
            ix=d['cell'][te]==c;error=float(abs(d['effect'][te][ix]-p['magnitude'][ix]).mean());be=float(abs(d['effect'][te][ix]-baseline[ix]).mean());skill=1-error/be
            metrics['cells'][str(d['cell_names'][c])].update(line_baseline_MAE=be,MAE_skill_vs_line_baseline=skill)
            if support[c]:skills.append(skill)
        metrics['macro_MAE_skill_vs_line_baseline_supported_cells']=float(np.mean(skills)) if skills else None
        curve.append({'epoch':epoch,'mean_data_objective':float(np.mean(objectives)),'mean_relative_reader_anchor':float(np.mean(anchors)),
           'elapsed_seconds':elapsed_before+time.monotonic()-start,**metrics});write_json(directory/(name+'_curve.json'),curve)
        # Save every prespecified epoch as a small override of the pinned base.
        state={k:v.detach().cpu().contiguous() for k,v in model.named_parameters() if k.startswith('head.') or (args.recipe=='adapt_last_stage' and k.startswith('reader.main.3.'))}
        checkpoint=directory/f'{name}_epoch{epoch}.safetensors';save_file(state,str(checkpoint));np.savez_compressed(directory/f'{name}_epoch{epoch}_predictions.npz',observation_index=te,**p)
        write_json(directory/f'{name}_epoch{epoch}_receipt.json',{'checkpoint_sha256':sha(checkpoint),'parameter_names':list(state),
           'prediction_sha256':sha(directory/f'{name}_epoch{epoch}_predictions.npz'),'common_epoch_selection_not_yet_applied':True})
        checkpoint_now('epoch_complete')
        print('NATIVE EPOCH COMPLETE',name,epoch,'magnitude_skill',metrics['macro_MAE_skill_vs_line_baseline_supported_cells'],flush=True)
    changes={k:int(torch.count_nonzero(v.detach().cpu()!=initial_reader[k])) for k,v in model.reader.state_dict().items()}
    assert sum(v for k,v in changes.items() if not k.startswith('main.3.'))==0
    assert sum(v for k,v in changes.items() if 'running_' in k or 'num_batches' in k)==0
    if args.recipe=='frozen_reader':assert sum(changes.values())==0
    else:assert sum(changes.values())>0
    write_json(record_path,{**meta,'completed_epochs':plan['epochs'],'actual_backprop_updates':updates,
       'actual_batch_schedule_sha256':schedule.hexdigest(),'reader_changed_scalar_counts':changes,
       'only_declared_reader_stage_can_change':True,'BN_running_statistics_unchanged':True,'scientific_goal_complete':False})
    assert_versions_stable();print('NATIVE FIT COMPLETE',name,flush=True)


if __name__=='__main__':main()

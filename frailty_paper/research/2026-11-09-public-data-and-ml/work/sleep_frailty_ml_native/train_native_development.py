"""Matched native-DNA effect-size fitting, with actual last-stage adaptation."""
import argparse,time
from native_common import *
from native_data import load_data,PackedDNA,normalization
from native_network import NativeAdaptationModel,torch
from training_utilities import transformations,detailed_metrics
from activity_effect_network import activity_loss
from safetensors.torch import save_file


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--recipe',required=True);parser.add_argument('--held-study',required=True);args=parser.parse_args()
    assert_versions_native();plan=get_plan();assert args.recipe in plan['recipes'] and args.held_study in plan['phase_one_roles']+plan['phase_two_roles_if_gate_passes']
    assert json.loads((DEST/'native_runtime_and_gradient_audit.json').read_text())['passed']
    d,pairs,bank,x,extra=load_data();added_siraj,added_neural=add_all_documented_signed_targets(d);packed=PackedDNA(d,bank);train=d['split']==0
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
    seed=plan['seed'];torch.manual_seed(seed);device='mps' if torch.backends.mps.is_available() else 'cpu'
    model=NativeAdaptationModel(med,lm,parent,norm,args.recipe,seed).to(device)
    model.anchored={k:v.to(device) for k,v in model.anchored.items()}
    initial_reader={k:v.detach().cpu().clone() for k,v in model.reader.state_dict().items()};h=hashlib.sha256()
    for key,value in sorted(model.state_dict().items()):h.update(key.encode());h.update(value.detach().cpu().numpy().tobytes())
    directory=DEST/'development';directory.mkdir(exist_ok=True);name=f'{args.recipe}_{args.held_study}_{seed}';record_path=directory/(name+'_record.json');assert not record_path.exists()
    assert not (directory/(name+'_started.json')).exists(),'Preserve and inspect any incomplete prior fit before a restart.'
    meta={'recipe':args.recipe,'held_study':args.held_study,'seed':seed,'epochs':plan['epochs'],'training_rows':len(tr),'prediction_rows':len(te),
       'training_index_sha256':hashlib.sha256(tr.tobytes()).hexdigest(),'prediction_index_sha256':hashlib.sha256(te.tobytes()).hexdigest(),
       'initial_parameter_sha256':h.hexdigest(),'normalization':norm,'cell_medians':med.tolist(),'line_medians':lm.tolist(),'line_parent':parent.tolist(),
       'supported_cells':support.tolist(),'supported_lines':ls.tolist(),'context_transformations':context,
       'new_original_Siraj_signed_training_rows':int(added_siraj[tr].sum()),'new_original_neuronal_signed_training_rows':int(added_neural[tr].sum()),
       'total_known_signed_training_rows':int(np.isfinite(d['signed_effect'][tr]).sum()),
       'trainable_head_parameters':sum(p.numel() for p in model.head.parameters()),
       'trainable_reader_parameters':sum(p.numel() for p in model.reader.parameters() if p.requires_grad),
       'plan_sha256':sha(DEST/'finite_native_adaptation_plan.json'),'pretrained_reader_sha256':plan['reader_sha256'],
       'code_sha256':{p.name:sha(p) for p in [Path(__file__),Path(__file__).with_name('native_network.py'),Path(__file__).with_name('native_data.py')]},
       'device':device,'torch_version':torch.__version__,'physical_source_receipt':packed.receipt,
       'source_standard_error_used_at_inference':False,'calibration_rows_used_for_fit':False,'reserved_Ulirsch_outcomes_used':False,'PID':os.getpid()}
    write_json(directory/(name+'_started.json'),meta);opt=plan['optimization']
    groups=[{'params':list(model.head.parameters()),'lr':opt['head_lr']}]
    core=[p for p in model.reader.parameters() if p.requires_grad]
    if core:groups.append({'params':core,'lr':opt['reader_lr']})
    optimizer=torch.optim.AdamW(groups,weight_decay=opt['weight_decay']);tp={k:p[np.all(train[p],axis=1)] for k,p in pairs.items()}
    baseline=np.where(ls[d['line'][te]],lm[d['line'][te]],med[d['cell'][te]])
    model.eval()
    with torch.no_grad():assert abs(model(te[:32],d,packed,bank,x,extra,support,ls)['magnitude'].cpu().numpy()-baseline[:32]).max()<1e-6
    def predict():
        model.eval();result={k:np.empty(len(te),dtype='float32') for k in ['magnitude','signed_effect','activity_reference','activity_alternate']}
        with torch.no_grad():
            for i in range(0,len(te),256):
                p=model(te[i:i+256],d,packed,bank,x,extra,support,ls)
                for k in result:result[k][i:i+256]=p[k].cpu().numpy()
        return result
    rng=np.random.default_rng(seed);schedule=hashlib.sha256();updates=0;curve=[];start=time.monotonic();last_print=start
    for epoch in range(1,plan['epochs']+1):
        torch.manual_seed(seed+epoch)
        if device=='mps':torch.mps.manual_seed(seed+epoch)
        model.train();order=rng.permutation(tr);objectives=[];anchors=[]
        for i in range(0,len(order),opt['main_batch']):
            main_ids=order[i:i+opt['main_batch']];chunks=[main_ids];counts=[]
            for key,n in [('local',opt['local_pairs']),('cross',opt['cross_pairs'])]:
                p=tp[key];count=n if len(p) else 0;counts.append(count)
                if count:chunks.append(p[rng.integers(len(p),size=count)].reshape(-1))
            ids=np.concatenate(chunks);schedule.update(ids.tobytes());optimizer.zero_grad();prediction=model(ids,d,packed,bank,x,extra,support,ls)
            batch={k:torch.from_numpy(d[k][ids]).to(device) for k in ['effect','signed_effect','activity_reference','activity_alternate']};batch['activity_center']=torch.from_numpy(centers[ids]).to(device)
            objective=activity_loss(prediction,batch,torch.from_numpy(w[ids]).to(device),torch.from_numpy(quality[ids]).to(device),len(main_ids),*counts)
            anchor=model.anchor_loss();loss=objective+anchor;assert torch.isfinite(loss);loss.backward()
            assert all(p.grad is None or torch.isfinite(p.grad).all() for p in model.parameters())
            torch.nn.utils.clip_grad_norm_([p for p in model.parameters() if p.requires_grad],opt['gradient_clip']);optimizer.step();updates+=1
            objectives.append(float(objective.detach().cpu()));anchors.append(float(anchor.detach().cpu()))
            if time.monotonic()-last_print>30:
                progress={'fit':name,'epoch':epoch,'epoch_rows_completed':min(i+len(main_ids),len(tr)),
                   'epoch_rows_total':len(tr),'updates':updates,'elapsed_seconds':time.monotonic()-start,'PID':os.getpid()}
                write_json(directory/(name+'_live_progress.json'),progress)
                print('NATIVE FIT PROGRESS',json.dumps(progress),flush=True);last_print=time.monotonic()
        p=predict();metrics=detailed_metrics(d,te,p,med,support,pairs);skills=[]
        for c in np.unique(d['cell'][te]):
            ix=d['cell'][te]==c;error=float(abs(d['effect'][te][ix]-p['magnitude'][ix]).mean());be=float(abs(d['effect'][te][ix]-baseline[ix]).mean());skill=1-error/be
            metrics['cells'][str(d['cell_names'][c])].update(line_baseline_MAE=be,MAE_skill_vs_line_baseline=skill)
            if support[c]:skills.append(skill)
        metrics['macro_MAE_skill_vs_line_baseline_supported_cells']=float(np.mean(skills)) if skills else None
        curve.append({'epoch':epoch,'mean_data_objective':float(np.mean(objectives)),'mean_relative_reader_anchor':float(np.mean(anchors)),
           'elapsed_seconds':time.monotonic()-start,**metrics});write_json(directory/(name+'_curve.json'),curve)
        # Save every prespecified epoch as a small override of the pinned base.
        state={k:v.detach().cpu().contiguous() for k,v in model.named_parameters() if k.startswith('head.') or (args.recipe=='adapt_last_stage' and k.startswith('reader.main.3.'))}
        checkpoint=directory/f'{name}_epoch{epoch}.safetensors';save_file(state,str(checkpoint));np.savez_compressed(directory/f'{name}_epoch{epoch}_predictions.npz',observation_index=te,**p)
        write_json(directory/f'{name}_epoch{epoch}_receipt.json',{'checkpoint_sha256':sha(checkpoint),'parameter_names':list(state),
           'prediction_sha256':sha(directory/f'{name}_epoch{epoch}_predictions.npz'),'common_epoch_selection_not_yet_applied':True})
        print('NATIVE EPOCH COMPLETE',name,epoch,'magnitude_skill',metrics['macro_MAE_skill_vs_line_baseline_supported_cells'],flush=True)
    changes={k:int(torch.count_nonzero(v.detach().cpu()!=initial_reader[k])) for k,v in model.reader.state_dict().items()}
    assert sum(v for k,v in changes.items() if not k.startswith('main.3.'))==0
    assert sum(v for k,v in changes.items() if 'running_' in k or 'num_batches' in k)==0
    if args.recipe=='frozen_reader':assert sum(changes.values())==0
    else:assert sum(changes.values())>0
    write_json(record_path,{**meta,'completed_epochs':plan['epochs'],'actual_backprop_updates':updates,
       'actual_batch_schedule_sha256':schedule.hexdigest(),'reader_changed_scalar_counts':changes,
       'only_declared_reader_stage_can_change':True,'BN_running_statistics_unchanged':True,'scientific_goal_complete':False})
    assert_versions_native();print('NATIVE FIT COMPLETE',name,flush=True)


if __name__=='__main__':main()

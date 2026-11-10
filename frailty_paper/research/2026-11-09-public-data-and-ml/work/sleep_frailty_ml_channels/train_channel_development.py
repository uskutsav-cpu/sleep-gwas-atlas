"""Matched finite comparison of summary-only and channel-preserving inputs."""
import argparse,time
from channel_common import *
from coverage_inputs import load_inputs
from training_utilities import transformations,detailed_metrics
from activity_effect_network import activity_loss
from cell_matched_effect_network import CellMatchedEffectHead,torch
from safetensors.torch import save_file

def channel_normalization(extra,banks):
    total=np.zeros(192,dtype='float64');squares=np.zeros(192,dtype='float64');n=0
    for i in range(0,len(banks),4096):
        values=extra[banks[i:i+4096]].reshape(-1,192).astype('float64');assert np.isfinite(values).all();total+=values.sum(0);squares+=(values*values).sum(0);n+=len(values)
    mean=total/n;sd=np.sqrt(np.maximum(squares/n-mean*mean,0));sd[sd<1e-5]=1
    return mean.astype('float32'),sd.astype('float32')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--recipe',choices=['scalar_mask','channel_preserving'],required=True);parser.add_argument('--held-study',required=True);args=parser.parse_args();assert_frozen_previous();plan=get_plan();assert args.held_study in plan['held_studies']
    receipt=json.loads((DEST/'channel_feature_receipt.json').read_text());assert receipt['passed'] and sha(DEST/'projected_reader_channels.npy')==receipt['features_sha256'];extra=np.load(DEST/'projected_reader_channels.npy',mmap_mode='r');d,pairs,bank,x=load_inputs();train=d['split']==0
    if args.held_study=='none':test=d['split']==1
    else:
        sid=list(d['study_names']).index(args.held_study);groups=np.unique(d['group'][d['study']==sid]);train&=(d['study']!=sid)&~np.isin(d['group'],groups);test=(d['study']==sid)&(d['split']!=2)
    tr,te=np.flatnonzero(train),np.flatnonzero(test);assert not set(d['group'][tr])&set(d['group'][te]);med,support,w,centers,quality,context=transformations(d,train,True)
    train_banks=np.unique(d['bank_index'][tr]);mean=x[train_banks].mean(0);sd=x[train_banks].std(0);sd[sd<1e-5]=1;emean,esd=channel_normalization(extra,train_banks)
    nc=len(med);nl=len(d['cell_line_names']);parent=np.zeros(nl,dtype='int64');line_med=np.zeros(nl,dtype='float32');line_support=np.zeros(nl,dtype=bool)
    for line in range(nl):
        cells=np.unique(d['cell'][d['line']==line]);assert len(cells)<=1
        if len(cells):parent[line]=cells[0]
        ix=train&(d['line']==line);line_support[line]=ix.any();line_med[line]=np.median(d['effect'][ix]) if ix.any() else med[parent[line]]
    seed=plan['seed'];torch.manual_seed(seed);model=CellMatchedEffectHead(med,line_med,parent,258,seed);h=hashlib.sha256()
    for key,value in sorted(model.state_dict().items()):h.update(key.encode());h.update(value.detach().numpy().tobytes())
    directory=DEST/'development';directory.mkdir(exist_ok=True);name=f'{args.recipe}_{args.held_study}_{seed}';record=directory/(name+'_record.json');assert not record.exists();metadata={'recipe':args.recipe,'held_study':args.held_study,'epochs':plan['epochs'],'seed':seed,'training_rows':len(tr),'prediction_rows':len(te),
        'training_index_sha256':hashlib.sha256(tr.tobytes()).hexdigest(),'prediction_index_sha256':hashlib.sha256(te.tobytes()).hexdigest(),'initial_parameter_sha256':h.hexdigest(),'cell_medians':med.tolist(),'line_medians':line_med.tolist(),'line_parent':parent.tolist(),'supported_cells':support.tolist(),'supported_lines':line_support.tolist(),
        'scalar_feature_mean':mean.tolist(),'scalar_feature_sd':sd.tolist(),'channel_feature_mean':emean.tolist(),'channel_feature_sd':esd.tolist(),'context_transformations':context,
        'trainable_head_parameters':sum(p.numel() for p in model.parameters()),'plan_sha256':sha(DEST/'finite_channel_plan.json'),'features_sha256':receipt['features_sha256'],'source_data_sha256':plan['source_data_sha256'],
        'pretrained_weights_updated':False,'previous_reserved_blood_or_OL33_data_used':False,'code_sha256':sha(Path(__file__)),'PID':os.getpid()}
    write_json(directory/(name+'_started.json'),metadata);optimizer=torch.optim.AdamW(model.parameters(),lr=3e-4,weight_decay=.001);tp={k:p[np.all(train[p],axis=1)] for k,p in pairs.items()}
    def forward(ids):
        b=d['bank_index'][ids];base=np.clip((x[b]-mean)/sd,-12,12);added=np.clip((extra[b].reshape(-1,192).astype('float32')-emean)/esd,-12,12)
        if args.recipe=='scalar_mask':added[:]=0
        features=np.concatenate([base,added],axis=1);raw=bank[b];c=d['cell'][ids];line=d['line'][ids];cf=np.eye(nc,dtype='float32')[c]*support[c,None];lf=np.eye(nl,dtype='float32')[line]*line_support[line,None]
        return model(torch.from_numpy(features),torch.from_numpy(raw[:,:,0].copy()),torch.from_numpy(raw[:,:,1].copy()),torch.from_numpy(c.astype('int64')),torch.from_numpy(cf),torch.from_numpy(lf))
    def predict():
        model.eval();output={k:np.empty(len(te),dtype='float32') for k in ['magnitude','signed_effect','activity_reference','activity_alternate']}
        with torch.no_grad():
            for i in range(0,len(te),2048):
                p=forward(te[i:i+2048])
                for key in output:output[key][i:i+2048]=p[key].numpy()
        return output
    baseline=np.where(line_support[d['line'][te]],line_med[d['line'][te]],med[d['cell'][te]]);model.eval()
    with torch.no_grad():assert abs(forward(te[:32])['magnitude'].numpy()-baseline[:32]).max()<1e-6
    curve=[];best=-np.inf;rng=np.random.default_rng(seed);batch_hash=hashlib.sha256();updates=0;start=time.monotonic()
    for epoch in range(1,plan['epochs']+1):
        torch.manual_seed(seed+epoch);model.train();order=rng.permutation(tr);losses=[]
        for i in range(0,len(order),1024):
            main_ids=order[i:i+1024];chunks=[main_ids];counts=[]
            for label in ['local','cross']:
                p=tp[label];count=32 if len(p) else 0;counts.append(count)
                if count:chunks.append(p[rng.integers(len(p),size=count)].reshape(-1))
            ids=np.concatenate(chunks);batch_hash.update(ids.tobytes());optimizer.zero_grad();prediction=forward(ids);batch={k:torch.from_numpy(d[k][ids]) for k in ['effect','signed_effect','activity_reference','activity_alternate']};batch['activity_center']=torch.from_numpy(centers[ids]);loss=activity_loss(prediction,batch,torch.from_numpy(w[ids]),torch.from_numpy(quality[ids]),len(main_ids),*counts);assert torch.isfinite(loss);loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),1.);optimizer.step();updates+=1;losses.append(float(loss.detach()))
        p=predict();metrics=detailed_metrics(d,te,p,med,support,pairs);skills=[]
        for c in np.unique(d['cell'][te]):
            ix=d['cell'][te]==c;error=float(abs(d['effect'][te][ix]-p['magnitude'][ix]).mean());be=float(abs(d['effect'][te][ix]-baseline[ix]).mean());skill=1-error/be;metrics['cells'][str(d['cell_names'][c])].update(line_baseline_MAE=be,MAE_skill_vs_line_baseline=skill)
            if support[c]:skills.append(skill)
        score=float(np.mean(skills)) if skills else -1e10;metrics['macro_MAE_skill_vs_line_baseline_supported_cells']=score;curve.append({'epoch':epoch,'mean_objective':float(np.mean(losses)),'elapsed_seconds':time.monotonic()-start,**metrics});write_json(directory/(name+'_curve.json'),curve)
        if score>best:
            best=score;checkpoint=directory/(name+'.safetensors');save_file({k:v.detach().contiguous() for k,v in model.state_dict().items()},checkpoint);np.savez_compressed(directory/(name+'_predictions.npz'),observation_index=te,**p);write_json(directory/(name+'_best.json'),{'epoch':epoch,'score':score,'checkpoint_sha256':sha(checkpoint)})
        print('CHANNEL DEVELOPMENT',name,epoch,'error_skill',round(score,6),flush=True)
    write_json(record,{**metadata,'completed_epochs':plan['epochs'],'actual_backprop_updates':updates,'actual_batch_schedule_sha256':batch_hash.hexdigest(),'checkpoint_sha256':sha(checkpoint),'best_development_skill':best,'scope':'Development only; no prior test is claimed untouched for this changed model.'});assert_frozen_previous();print('CHANNEL FIT COMPLETE',name,flush=True)

if __name__=='__main__':main()

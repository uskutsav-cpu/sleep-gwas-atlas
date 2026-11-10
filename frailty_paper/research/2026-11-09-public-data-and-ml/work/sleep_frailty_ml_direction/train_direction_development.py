"""Fit actual observed magnitudes with and without declared measurement auxiliaries."""
import argparse,time
from direction_common import *
from measurement_inputs import load_data,normalization,normalized_inputs
from finite_direction_plan import get_plan
from folded_effect_network import FoldedEffectHead,measurement_auxiliary,torch
from activity_effect_network import activity_loss
from training_utilities import transformations,detailed_metrics
from safetensors.torch import save_file

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--recipe',required=True);parser.add_argument('--held-study',required=True);args=parser.parse_args();assert_versions();plan=get_plan();assert args.recipe in plan['recipes'] and args.held_study in plan['held_studies'];assert json.loads((MEASUREMENT/'folded_median_pretraining_audit.json').read_text())['passed']
    d,pairs,bank,x,extra=load_data();new_signed=load_signed_sidecar(d);train=d['split']==0
    if args.held_study=='none':test=d['split']==1
    else:
        sid=list(d['study_names']).index(args.held_study);groups=np.unique(d['group'][d['study']==sid]);train&=(d['study']!=sid)&~np.isin(d['group'],groups);test=(d['study']==sid)&(d['split']!=2)
    tr,te=np.flatnonzero(train),np.flatnonzero(test);assert not set(d['group'][tr])&set(d['group'][te]);med,support,weights,centers,quality,context=transformations(d,train,True);norm=normalization(extra,x,np.unique(d['bank_index'][tr]));nc=len(med);nl=len(d['cell_line_names']);parent=np.zeros(nl,dtype='int64');line_med=np.zeros(nl,dtype='float32');line_support=np.zeros(nl,dtype=bool)
    for line in range(nl):
        cells=np.unique(d['cell'][d['line']==line]);assert len(cells)<=1
        if len(cells):parent[line]=cells[0]
        ix=train&(d['line']==line);line_support[line]=ix.any();line_med[line]=np.median(d['effect'][ix]) if ix.any() else med[parent[line]]
    seed=plan['seed'];torch.manual_seed(seed);model=FoldedEffectHead(med,line_med,parent,seed);h=hashlib.sha256()
    for key,value in sorted(model.state_dict().items()):h.update(key.encode());h.update(value.detach().numpy().tobytes())
    out=DEST/'development';out.mkdir(exist_ok=True);name=f'{args.recipe}_{args.held_study}_{seed}';assert not (out/(name+'_record.json')).exists();metadata={'recipe':args.recipe,'held_study':args.held_study,'seed':seed,'epochs':plan['epochs'],'training_rows':len(tr),'prediction_rows':len(te),
        'training_index_sha256':hashlib.sha256(tr.tobytes()).hexdigest(),'prediction_index_sha256':hashlib.sha256(te.tobytes()).hexdigest(),'initial_parameter_sha256':h.hexdigest(),'cell_medians':med.tolist(),'line_medians':line_med.tolist(),'line_parent':parent.tolist(),'supported_cells':support.tolist(),'supported_lines':line_support.tolist(),
        'normalization':norm,'context_transformations':context,'trainable_parameters':sum(p.numel() for p in model.parameters()),'source_data_sha256':plan['data_sha256'],'plan_sha256':sha(DEST/'finite_direction_plan.json'),
        'code_sha256':{p.name:sha(p) for p in [Path(__file__),MEASURE_WORK/'folded_effect_network.py',MEASURE_WORK/'measurement_inputs.py']},'new_original_signed_training_targets':int(new_signed[tr].sum()),'new_Ulirsch_outcomes_used':False,'pretrained_weights_updated':False,'PID':os.getpid()}
    write_json(out/(name+'_started.json'),metadata);optimizer=torch.optim.AdamW(model.parameters(),lr=3e-4,weight_decay=.001);tp={k:p[np.all(train[p],axis=1)] for k,p in pairs.items()}
    def forward(ids):
        b=d['bank_index'][ids];mean,delta,physical=normalized_inputs(extra,x,b,norm);raw=bank[b];c=d['cell'][ids];line=d['line'][ids];cf=np.eye(nc,dtype='float32')[c]*support[c,None];lf=np.eye(nl,dtype='float32')[line]*line_support[line,None]
        return model(torch.from_numpy(mean),torch.from_numpy(delta),torch.from_numpy(physical),torch.from_numpy(raw[:,:,0].copy()),torch.from_numpy(raw[:,:,1].copy()),torch.from_numpy(c.astype('int64')),torch.from_numpy(cf),torch.from_numpy(lf))
    def predict():
        model.eval();result={k:np.empty(len(te),dtype='float32') for k in ['magnitude','signed_effect','activity_reference','activity_alternate','assay_SD_component','latent_absolute_component']}
        with torch.no_grad():
            for i in range(0,len(te),2048):
                p=forward(te[i:i+2048])
                for k in result:result[k][i:i+2048]=p[k].numpy()
        return result
    baseline=np.where(line_support[d['line'][te]],line_med[d['line'][te]],med[d['cell'][te]]);curve=[];best=-np.inf;rng=np.random.default_rng(seed);schedule=hashlib.sha256();updates=0;start=time.monotonic()
    for epoch in range(1,plan['epochs']+1):
        torch.manual_seed(seed+epoch);model.train();order=rng.permutation(tr);losses=[];main_losses=[];extras=[]
        for i in range(0,len(order),1024):
            main=order[i:i+1024];chunks=[main];counts=[]
            for key in ['local','cross']:
                p=tp[key];count=32 if len(p) else 0;counts.append(count)
                if count:chunks.append(p[rng.integers(len(p),size=count)].reshape(-1))
            ids=np.concatenate(chunks);schedule.update(ids.tobytes());optimizer.zero_grad();p=forward(ids);batch={k:torch.from_numpy(d[k][ids]) for k in ['effect','signed_effect','activity_reference','activity_alternate']};batch['activity_center']=torch.from_numpy(centers[ids]);loss=activity_loss(p,batch,torch.from_numpy(weights[ids]),torch.from_numpy(quality[ids]),len(main),*counts);main_losses.append(float(loss.detach()));aux=loss*0
            if args.recipe=='measurement_consistency':aux=measurement_auxiliary(p,batch,torch.from_numpy(d['measurement_standard_error'][ids]),torch.from_numpy(weights[ids]),torch.from_numpy(quality[ids]),len(main))
            loss=loss+aux;assert torch.isfinite(loss);loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),1.);optimizer.step();updates+=1;losses.append(float(loss.detach()));extras.append(float(aux.detach()))
        p=predict();metrics=detailed_metrics(d,te,p,med,support,pairs);skills=[]
        for c in np.unique(d['cell'][te]):
            ix=d['cell'][te]==c;error=float(abs(d['effect'][te][ix]-p['magnitude'][ix]).mean());be=float(abs(d['effect'][te][ix]-baseline[ix]).mean());skill=1-error/be;metrics['cells'][str(d['cell_names'][c])].update(line_baseline_MAE=be,MAE_skill_vs_line_baseline=skill)
            if support[c]:skills.append(skill)
        score=float(np.mean(skills)) if skills else -1e10;metrics['macro_MAE_skill_vs_line_baseline_supported_cells']=score;curve.append({'epoch':epoch,'mean_total_objective':float(np.mean(losses)),'mean_primary_RNA_pair_objective':float(np.mean(main_losses)),'mean_measurement_auxiliary':float(np.mean(extras)),'elapsed_seconds':time.monotonic()-start,**metrics});write_json(out/(name+'_curve.json'),curve)
        if score>best:
            best=score;checkpoint=out/(name+'.safetensors');save_file({k:v.detach().contiguous() for k,v in model.state_dict().items()},checkpoint);np.savez_compressed(out/(name+'_predictions.npz'),observation_index=te,**p);write_json(out/(name+'_best.json'),{'epoch':epoch,'score':score,'checkpoint_sha256':sha(checkpoint)})
        print('MEASUREMENT DEVELOPMENT',name,epoch,'observed_magnitude_error_skill',round(score,6),flush=True)
    write_json(out/(name+'_record.json'),{**metadata,'completed_epochs':plan['epochs'],'actual_backprop_updates':updates,'actual_batch_schedule_sha256':schedule.hexdigest(),'checkpoint_sha256':sha(checkpoint),'best_development_skill':best,'source_SD_used_at_inference':False,'latent_or_SD_component_not_claimed_validated':True});assert_versions();print('MEASUREMENT FIT COMPLETE',name,flush=True)

if __name__=='__main__':main()

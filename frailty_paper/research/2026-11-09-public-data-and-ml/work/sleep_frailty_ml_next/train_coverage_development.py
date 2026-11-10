"""Finite single-architecture coverage expansion; no reserved predictions or labels."""
import sys,argparse,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
from training_utilities import transformations,detailed_metrics
from activity_effect_network import activity_loss
from cell_matched_effect_network import CellMatchedEffectHead,torch
from coverage_inputs import load_inputs,DEST
from safetensors.torch import save_file

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--held-study',required=True);args=parser.parse_args()
    assert_previous_sources();assert json.loads((DEST/'coverage_independent_data_audit.json').read_text())['passed']
    plan=json.loads((DEST/'coverage_development_plan.json').read_text());assert args.held_study in plan['held_studies'];seed=plan['seed'];epochs=plan['epochs']
    d,pairs,bank,x=load_inputs();train=d['split']==0
    if args.held_study=='none':test=d['split']==1
    else:
        sid=list(d['study_names']).index(args.held_study);groups=np.unique(d['group'][d['study']==sid]);train&=(d['study']!=sid)&~np.isin(d['group'],groups);test=(d['study']==sid)&(d['split']!=2)
    tr,te=np.flatnonzero(train),np.flatnonzero(test);assert not set(d['group'][tr])&set(d['group'][te])
    med,support,w,centers,quality,context=transformations(d,train,True);train_banks=np.unique(d['bank_index'][tr]);mean=x[train_banks].mean(0);sd=x[train_banks].std(0);sd[sd<1e-5]=1
    nl=len(d['cell_line_names']);nc=len(med);parent=np.zeros(nl,dtype='int64');line_med=np.zeros(nl,dtype='float32');line_support=np.zeros(nl,dtype=bool)
    for line in range(nl):
        cells=np.unique(d['cell'][d['line']==line]);assert len(cells)<=1
        if len(cells):parent[line]=cells[0]
        ix=train&(d['line']==line);line_support[line]=ix.any();line_med[line]=np.median(d['effect'][ix]) if ix.any() else med[parent[line]]
    torch.manual_seed(seed);model=CellMatchedEffectHead(med,line_med,parent,66,seed);h=hashlib.sha256()
    for k,v in sorted(model.state_dict().items()):h.update(k.encode());h.update(v.detach().numpy().tobytes())
    optimizer=torch.optim.AdamW(model.parameters(),lr=3e-4,weight_decay=.001);output=DEST/'coverage_development';output.mkdir(exist_ok=True);name=f'three_readers_{args.held_study}_{seed}';assert not (output/(name+'_record.json')).exists()
    metadata={'recipe':'three_readers','held_study':args.held_study,'epochs':epochs,'seed':seed,'training_rows':len(tr),'prediction_rows':len(te),
              'training_index_sha256':hashlib.sha256(tr.tobytes()).hexdigest(),'prediction_index_sha256':hashlib.sha256(te.tobytes()).hexdigest(),'initial_head_sha256':h.hexdigest(),
              'cell_medians':med.tolist(),'line_medians':line_med.tolist(),'line_parent':parent.tolist(),'supported_cells':support.tolist(),'supported_lines':line_support.tolist(),
              'feature_mean':mean.tolist(),'feature_sd':sd.tolist(),'context_transformations':context,'trainable_head_parameters':sum(p.numel() for p in model.parameters()),
              'pretrained_reader_parameters_updated':False,'plan_sha256':sha(DEST/'coverage_development_plan.json'),'source_training_arrays_sha256':sha(DEST/'coverage_training_arrays.npz'),
              'code_sha256':{p.name:sha(p) for p in [Path(__file__),Path(__file__).with_name('cell_matched_effect_network.py'),Path(__file__).with_name('coverage_inputs.py')]},
              'reserved_inputs_or_outcomes_used':False,'PID':os.getpid()}
    write_json(output/(name+'_started.json'),metadata);tp={k:p[np.all(train[p],axis=1)] for k,p in pairs.items()}
    def forward(ids):
        b=d['bank_index'][ids];features=np.clip((x[b]-mean)/sd,-12,12);raw=bank[b];c=d['cell'][ids];line=d['line'][ids]
        cf=np.eye(nc,dtype='float32')[c]*support[c,None];lf=np.eye(nl,dtype='float32')[line]*line_support[line,None]
        return model(torch.from_numpy(features),torch.from_numpy(raw[:,:,0].copy()),torch.from_numpy(raw[:,:,1].copy()),torch.from_numpy(c.astype('int64')),torch.from_numpy(cf),torch.from_numpy(lf))
    def predict():
        model.eval();p={k:np.empty(len(te),dtype='float32') for k in ['magnitude','signed_effect','activity_reference','activity_alternate']}
        with torch.no_grad():
            for i in range(0,len(te),2048):
                got=forward(te[i:i+2048])
                for k in p:p[k][i:i+2048]=got[k].numpy()
        return p
    baseline=np.where(line_support[d['line'][te]],line_med[d['line'][te]],med[d['cell'][te]]);model.eval()
    with torch.no_grad():initial=forward(te[:32])['magnitude'].numpy()
    assert abs(initial-baseline[:32]).max()<1e-6
    curve=[];best=-np.inf;updates=0;rng=np.random.default_rng(seed);started=time.monotonic()
    for epoch in range(1,epochs+1):
        torch.manual_seed(seed+epoch);model.train();order=rng.permutation(tr);losses=[]
        for i in range(0,len(order),1024):
            main_ids=order[i:i+1024];chunks=[main_ids];counts=[]
            for key in ['local','cross']:
                p=tp[key];count=32 if len(p) else 0;counts.append(count)
                if count:chunks.append(p[rng.integers(len(p),size=count)].reshape(-1))
            ids=np.concatenate(chunks);optimizer.zero_grad();prediction=forward(ids);batch={k:torch.from_numpy(d[k][ids]) for k in ['effect','signed_effect','activity_reference','activity_alternate']};batch['activity_center']=torch.from_numpy(centers[ids]);loss=activity_loss(prediction,batch,torch.from_numpy(w[ids]),torch.from_numpy(quality[ids]),len(main_ids),*counts)
            assert torch.isfinite(loss);loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),1.);optimizer.step();updates+=1;losses.append(float(loss.detach()))
        prediction=predict();metrics=detailed_metrics(d,te,prediction,med,support,pairs);skills=[]
        for c in np.unique(d['cell'][te]):
            ix=d['cell'][te]==c;error=float(abs(d['effect'][te][ix]-prediction['magnitude'][ix]).mean());base=float(abs(d['effect'][te][ix]-baseline[ix]).mean());skill=1-error/base;metrics['cells'][str(d['cell_names'][c])].update(line_baseline_MAE=base,MAE_skill_vs_line_baseline=skill)
            if support[c]:skills.append(skill)
        score=float(np.mean(skills)) if skills else -1e10;metrics['macro_MAE_skill_vs_line_baseline_supported_cells']=score;curve.append({'epoch':epoch,'mean_objective':float(np.mean(losses)),'elapsed_seconds':time.monotonic()-started,**metrics});write_json(output/(name+'_curve.json'),curve)
        if score>best:
            best=score;checkpoint=output/(name+'.safetensors');save_file({k:v.detach().contiguous() for k,v in model.state_dict().items()},checkpoint);np.savez_compressed(output/(name+'_predictions.npz'),observation_index=te,**prediction);write_json(output/(name+'_best.json'),{'epoch':epoch,'score':score,'checkpoint_sha256':sha(checkpoint)})
        print('COVERAGE DEVELOPMENT',name,epoch,'error_skill_vs_line_median',round(score,5),flush=True)
    write_json(output/(name+'_record.json'),{**metadata,'completed_epochs':epochs,'actual_backprop_updates':updates,'best_development_line_baseline_skill':best,'checkpoint_sha256':sha(checkpoint),
              'scope':'Selected development only. No independent test or disease-score validation established.'});print('COVERAGE FIT COMPLETE',name,flush=True)

if __name__=='__main__':main()

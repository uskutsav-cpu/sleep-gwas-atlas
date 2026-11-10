"""One declared line-aware development architecture, frozen versus tuned."""
import argparse,time
from v9_common import *
from line_aware_network import LineAwareEffectNetwork,torch,F
from feature_banks import global_bank,local_bank
from training_utilities import transformations,detailed_metrics
from activity_effect_network import activity_loss
from safetensors.torch import save_file


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--recipe',choices=['frozen','tuned'],required=True);parser.add_argument('--held-study',default='none');parser.add_argument('--epochs',type=int,default=8);parser.add_argument('--seed',type=int,default=42);args=parser.parse_args()
    assert_previous_sources();assert json.loads((OUT/'independent_augmented_partition_audit.json').read_text())['passed']
    with np.load(DATA/'training_arrays.npz') as z:d={k:z[k] for k in z.files}
    with np.load(DATA/'measured_pairs.npz') as z:pairs={k:z[k] for k in z.files}
    train_mask=d['split']==0
    if args.held_study=='none':test_mask=d['split']==1
    else:
        s=list(d['study_names']).index(args.held_study);groups=np.unique(d['group'][d['study']==s]);train_mask&=(d['study']!=s)&~np.isin(d['group'],groups);test_mask=(d['study']==s)&(d['split']!=2)
    tr,te=np.flatnonzero(train_mask),np.flatnonzero(test_mask);assert not set(d['group'][tr])&set(d['group'][te])
    med,support,weights,centers,quality,metadata=transformations(d,train_mask,True)
    nl=len(d['cell_line_names']);line_parent=np.zeros(nl,dtype='int64');line_med=np.full(nl,np.median(d['effect'][tr]),dtype='float32');line_support=np.zeros(nl,dtype=bool)
    for line in range(nl):
        allrows=d['line']==line;cells=np.unique(d['cell'][allrows])
        if len(cells):assert len(cells)==1;line_parent[line]=cells[0]
        ids=train_mask&allrows
        if ids.any():line_med[line]=np.median(d['effect'][ids]);line_support[line]=True
        else:line_med[line]=med[line_parent[line]]
    h,local=global_bank(),local_bank();torch.manual_seed(args.seed);model=LineAwareEffectNetwork(med,line_med,line_parent,args.recipe,args.seed)
    optimizer=torch.optim.AdamW([p for p in model.parameters() if p.requires_grad],lr=3e-4,weight_decay=.001)
    name=f'{args.recipe}_{args.held_study}_{args.seed}';assert not (OUT/(name+'_record.json')).exists();checkpoint=MODELS/(name+'.safetensors')
    info={'recipe':args.recipe,'held_study':args.held_study,'seed':args.seed,'epochs':args.epochs,'training_rows':len(tr),'test_rows':len(te),'data_sha256':sha(DATA/'training_arrays.npz'),'cell_medians':med.tolist(),'line_medians':line_med.tolist(),'line_parent_cell':line_parent.tolist(),'supported_cells':support.tolist(),'supported_lines':line_support.tolist(),'activity_and_precision_transformations':metadata,'training_index_sha256':hashlib.sha256(tr.tobytes()).hexdigest(),'test_index_sha256':hashlib.sha256(te.tobytes()).hexdigest(),'PID':os.getpid(),'fresh_Matoba_model_performance_read':False,'code_sha256':{p.name:sha(p) for p in [Path(__file__),WORK/'line_aware_network.py']}}
    write_json(OUT/(name+'_started.json'),info)
    if not (OUT/'model_search_started.json').exists():write_json(OUT/'model_search_started.json',{'amendment_plan_sha256':sha(OUT/'revealed_source_amendment_plan.json'),'fresh_test_still_unrevealed':True,'PID':os.getpid()})
    tp={k:p[np.all(train_mask[p],axis=1)] for k,p in pairs.items()}
    def forward(ids):
        banks=d['bank_index'][ids];gg=torch.from_numpy(h[banks]);ll=torch.from_numpy(local[banks]);c=d['cell'][ids];line=d['line'][ids]
        cf=np.eye(len(med),dtype='float32')[c];cf*=support[c,None];lf=np.eye(nl,dtype='float32')[line];lf*=line_support[line,None]
        return model(gg[:,0],gg[:,1],ll[:,0],ll[:,1],torch.from_numpy(c.astype('int64')),torch.from_numpy(line.astype('int64')),torch.from_numpy(cf),torch.from_numpy(lf))
    def predict():
        model.eval();out={k:np.empty(len(te),dtype='float32') for k in ['magnitude','signed_effect','activity_reference','activity_alternate']}
        with torch.no_grad():
            for i in range(0,len(te),2048):
                got=forward(te[i:i+2048])
                for k in out:out[k][i:i+2048]=got[k].numpy()
        return out
    baseline=np.where(line_support[d['line'][te]],line_med[d['line'][te]],med[d['cell'][te]])
    model.eval()
    with torch.no_grad():
        initial=forward(te[:16])['magnitude'].numpy()
    assert np.max(abs(initial-baseline[:16]))<1e-6
    rng=np.random.default_rng(args.seed);curve=[];best=-np.inf;updates=0;start=time.monotonic()
    for epoch in range(1,args.epochs+1):
        torch.manual_seed(args.seed+epoch);model.train();order=rng.permutation(tr);losses=[]
        for i in range(0,len(order),1024):
            main_ids=order[i:i+1024];chunks=[main_ids];counts=[]
            for key in ['local','cross']:
                p=tp[key];count=32 if len(p) else 0;counts.append(count)
                if count:chunks.append(p[rng.integers(len(p),size=count)].reshape(-1))
            ids=np.concatenate(chunks);optimizer.zero_grad();pred=forward(ids)
            batch={k:torch.from_numpy(d[k][ids]) for k in ['effect','signed_effect','activity_reference','activity_alternate']};batch['activity_center']=torch.from_numpy(centers[ids])
            loss=activity_loss(pred,batch,torch.from_numpy(weights[ids]),torch.from_numpy(quality[ids]),len(main_ids),*counts);assert torch.isfinite(loss);loss.backward();torch.nn.utils.clip_grad_norm_([p for p in model.parameters() if p.requires_grad],1.);optimizer.step();updates+=1;losses.append(float(loss.detach()))
        pred=predict();met=detailed_metrics(d,te,pred,med,support,pairs);skills=[]
        for c in np.unique(d['cell'][te]):
            ix=d['cell'][te]==c;error=float(abs(pred['magnitude'][ix]-d['effect'][te][ix]).mean());base=float(abs(baseline[ix]-d['effect'][te][ix]).mean());skill=1-error/base
            met['cells'][str(d['cell_names'][c])].update(line_baseline_MAE=base,MAE_skill_vs_line_baseline=skill,same_line_supported_fraction=float(line_support[d['line'][te][ix]].mean()))
            if support[c]:skills.append(skill)
        score=float(np.mean(skills)) if skills else -1e10
        met['macro_MAE_skill_vs_line_baseline_supported_cells']=score
        curve.append({'epoch':epoch,'training_objective':float(np.mean(losses)),'elapsed_seconds':time.monotonic()-start,**met});write_json(OUT/(name+'_curve.json'),curve)
        if score>best:
            best=score;save_file({k:v.detach().contiguous() for k,v in model.state_dict().items()},checkpoint);write_json(OUT/(name+'_best_development.json'),{'epoch':epoch,'metrics':met,'checkpoint_sha256':sha(checkpoint)});np.savez_compressed(OUT/(name+'_predictions.npz'),observation_index=te,**pred)
        print('LINE-AWARE DEVELOPMENT',name,epoch,'error_skill_over_line_median',round(score,5),'seconds',round(time.monotonic()-start,1),flush=True)
    write_json(OUT/(name+'_record.json'),{**info,'completed_epochs':args.epochs,'actual_backprop_updates':updates,'checkpoint':str(checkpoint),'checkpoint_sha256':sha(checkpoint),'best_line_baseline_error_skill':best,'earlier_external_sources_are_revealed_development':True,'scope':'Development only. Fresh independent molecular and matched AVI tests remain; no genetic-score or disease validation.'})
    assert_previous_sources();print('LINE-AWARE DEVELOPMENT COMPLETE',name,flush=True)


if __name__=='__main__':main()

"""Fixed-epoch fits with matched shared parameters and withheld source predictions."""
import sys,argparse,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
from training_utilities import transformations
from activity_effect_network import activity_loss
from coverage_inputs import load_inputs,DEST
from matched_coverage_head import MatchedCoverageHead,torch
from safetensors.torch import save_file
FINAL=DEST/'validation'

def shared_hash(model,avi):
    h=hashlib.sha256()
    for key,value in sorted(model.state_dict().items()):
        if avi and key=='hidden.weight':value=torch.cat([value[:,:66],value[:,69:]],1)
        a=value.detach().numpy().copy();h.update(key.encode());h.update(str(a.shape).encode());h.update(a.tobytes())
    return h.hexdigest()

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--held-study',required=True);parser.add_argument('--seed',required=True,type=int);parser.add_argument('--avi',action='store_true');args=parser.parse_args()
    assert_previous_sources();lock_path=FINAL/'final_settings_lock.json';lock=json.loads(lock_path.read_text());assert args.seed in lock['seeds'] and args.held_study in lock['held_studies']
    for name,h in lock['fit_code_sha256'].items():assert sha(Path(__file__).parent/name)==h
    d,pairs,bank,x=load_inputs();assert sha(DEST/'coverage_training_arrays.npz')==lock['source_data_sha256'];train=d['split']!=2;other=None
    if args.held_study=='none':test=d['split']==2
    else:
        sid=list(d['study_names']).index(args.held_study);groups=np.unique(d['group'][d['study']==sid]);train&=(d['study']!=sid)&~np.isin(d['group'],groups);test=d['study']==sid
        other=(d['split']==2)&(d['study']!=sid)&~np.isin(d['group'],groups)
    tr,te=np.flatnonzero(train),np.flatnonzero(test);assert not set(d['group'][tr])&set(d['group'][te]);med,support,w,centers,quality,context=transformations(d,train,True)
    train_banks=np.unique(d['bank_index'][tr]);mean=x[train_banks].mean(0);sd=x[train_banks].std(0);sd[sd<1e-5]=1;nc=len(med);nl=len(d['cell_line_names']);parent=np.zeros(nl,dtype='int64');line_med=np.zeros(nl,dtype='float32');line_support=np.zeros(nl,dtype=bool)
    cell_sources=np.zeros(nc,dtype='int16');line_sources=np.zeros(nl,dtype='int16')
    for c in range(nc):cell_sources[c]=len(np.unique(d['study'][train&(d['cell']==c)]))
    for line in range(nl):
        cells=np.unique(d['cell'][d['line']==line]);assert len(cells)<=1
        if len(cells):parent[line]=cells[0]
        ix=train&(d['line']==line);line_support[line]=ix.any();line_med[line]=np.median(d['effect'][ix]) if ix.any() else med[parent[line]];line_sources[line]=len(np.unique(d['study'][ix]))
    av=None;transform=None
    if args.avi:
        receipt=json.loads((DEST/'coverage_AVI_input_receipt.json').read_text());assert receipt['passed'];assert sha(DEST/'coverage_appended_AVI.npy')==receipt['appended_bank_sha256']
        from feature_banks import avi_bank
        old=avi_bank();av=np.concatenate([old[np.arange(old.shape[0])],np.load(DEST/'coverage_appended_AVI.npy')]);values=av[train_banks];known=values[:,2]==1;assert known.any()
        median=np.nanmedian(values[known,:2],axis=0);filled=np.where(np.isfinite(values[:,:2]),values[:,:2],median);amean,asd=filled.mean(0),filled.std(0);asd[asd<1e-5]=1
        transform={'median':median.tolist(),'mean':amean.tolist(),'sd':asd.tolist()}
    torch.manual_seed(args.seed);model=MatchedCoverageHead(med,line_med,parent,args.seed,args.avi);arm='AVI' if args.avi else 'sequence';name=f'fixed_{args.held_study}_{arm}_{args.seed}';directory=FINAL/'fits';directory.mkdir(exist_ok=True);record=directory/(name+'_record.json');assert not record.exists()
    checkpoint=directory/(name+'.safetensors');metadata={'arm':arm,'seed':args.seed,'held_study':args.held_study,'epochs':lock['epochs'],'training_rows':len(tr),'prediction_rows':len(te),
              'cell_medians':med.tolist(),'line_medians':line_med.tolist(),'line_parent':parent.tolist(),'supported_cells':support.tolist(),'supported_lines':line_support.tolist(),
              'cell_names':d['cell_names'].tolist(),'cell_line_names':d['cell_line_names'].tolist(),'training_cell_source_counts':cell_sources.tolist(),'training_line_source_counts':line_sources.tolist(),
              'feature_mean':mean.tolist(),'feature_sd':sd.tolist(),'activity_and_precision_transformations':context,'AVI_transform':transform,
              'shared_initial_parameter_sha256':shared_hash(model,args.avi),'training_index_sha256':hashlib.sha256(tr.tobytes()).hexdigest(),'prediction_index_sha256':hashlib.sha256(te.tobytes()).hexdigest(),
              'settings_lock_sha256':sha(lock_path),'code_sha256':sha(Path(__file__)),'trainable_head_parameters':sum(p.numel() for p in model.parameters()),
              'pretrained_reader_parameters_updated':False,'reserved_inputs_or_outcomes_used':False,'PID':os.getpid()}
    write_json(directory/(name+'_started.json'),metadata);optimizer=torch.optim.AdamW(model.parameters(),lr=3e-4,weight_decay=.001);tp={k:p[np.all(train[p],axis=1)] for k,p in pairs.items()}
    def forward(ids):
        b=d['bank_index'][ids];features=np.clip((x[b]-mean)/sd,-12,12);raw=bank[b];c=d['cell'][ids];line=d['line'][ids];cf=np.eye(nc,dtype='float32')[c]*support[c,None];lf=np.eye(nl,dtype='float32')[line]*line_support[line,None]
        extra=None
        if args.avi:
            a=av[b].copy();a[:,:2]=(np.where(np.isfinite(a[:,:2]),a[:,:2],median)-amean)/asd;extra=torch.from_numpy(np.clip(a,-12,12))
        return model(torch.from_numpy(features),torch.from_numpy(raw[:,:,0].copy()),torch.from_numpy(raw[:,:,1].copy()),torch.from_numpy(c.astype('int64')),torch.from_numpy(cf),torch.from_numpy(lf),extra)
    rng=np.random.default_rng(args.seed);updates=0;start=time.monotonic();batch_hash=hashlib.sha256()
    for epoch in range(1,lock['epochs']+1):
        torch.manual_seed(args.seed+epoch);model.train();order=rng.permutation(tr);losses=[]
        for i in range(0,len(order),1024):
            main_ids=order[i:i+1024];chunks=[main_ids];counts=[]
            for key in ['local','cross']:
                p=tp[key];count=32 if len(p) else 0;counts.append(count)
                if count:chunks.append(p[rng.integers(len(p),size=count)].reshape(-1))
            ids=np.concatenate(chunks);batch_hash.update(ids.tobytes());optimizer.zero_grad();prediction=forward(ids);batch={k:torch.from_numpy(d[k][ids]) for k in ['effect','signed_effect','activity_reference','activity_alternate']};batch['activity_center']=torch.from_numpy(centers[ids]);loss=activity_loss(prediction,batch,torch.from_numpy(w[ids]),torch.from_numpy(quality[ids]),len(main_ids),*counts)
            assert torch.isfinite(loss);loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),1.);optimizer.step();updates+=1;losses.append(float(loss.detach()))
        print('FIXED COVERAGE',name,epoch,round(time.monotonic()-start,1),flush=True)
    model.eval()
    def predict(indices):
        output={k:np.empty(len(indices),dtype='float32') for k in ['magnitude','signed_effect','activity_reference','activity_alternate']}
        with torch.no_grad():
            for i in range(0,len(indices),2048):
                got=forward(indices[i:i+2048])
                for k in output:output[k][i:i+2048]=got[k].numpy()
        return {'observation_index':indices,'original_split':d['split'][indices],**output}
    save_file({k:v.detach().contiguous() for k,v in model.state_dict().items()},checkpoint);predictions=directory/(name+'_predictions.npz');np.savez_compressed(predictions,**predict(te));extra_record={}
    if other is not None:
        oi=np.flatnonzero(other);assert not set(d['group'][tr])&set(d['group'][oi]);path=directory/(name+'_other_calibration_predictions.npz');np.savez_compressed(path,**predict(oi));extra_record={'other_calibration_predictions':str(path),'other_calibration_predictions_sha256':sha(path),'other_calibration_rows':len(oi)}
    write_json(record,{**metadata,**extra_record,'actual_backprop_updates':updates,'actual_batch_schedule_sha256':batch_hash.hexdigest(),'checkpoint':str(checkpoint),'checkpoint_sha256':sha(checkpoint),
              'predictions':str(predictions),'predictions_sha256':sha(predictions),'completed_epochs':lock['epochs'],'no_calibration_or_reserved_checkpoint_selection':True,'elapsed_seconds':time.monotonic()-start});assert_previous_sources();print('FIXED COVERAGE COMPLETE',name,flush=True)

if __name__=='__main__':main()

"""Matched fixed V9 fits, with no checkpoint selection on calibration or test."""
import argparse,time
from v9_common import *
from line_aware_network import LineAwareEffectNetwork,torch
from feature_banks import global_bank,local_bank,avi_bank
from training_utilities import transformations
from activity_effect_network import activity_loss
from safetensors.torch import save_file

FINAL=ROOT/'outputs/sleep_frailty_ai_v9_validation'


def shared_initial_hash(model,avi):
    h=hashlib.sha256()
    for k,v in sorted(model.state_dict().items()):
        if avi and k=='mag_hidden.weight':v=v[:,:-3]
        a=v.detach().numpy().copy();h.update(k.encode());h.update(str(a.shape).encode());h.update(a.tobytes())
    return h.hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--held-study',required=True);parser.add_argument('--seed',type=int,required=True);parser.add_argument('--avi',action='store_true');args=parser.parse_args()
    assert_previous_sources();lock_path=FINAL/'final_settings_lock.json';lock=json.loads(lock_path.read_text())
    assert args.seed in lock['seeds'] and args.held_study in lock['held_studies']
    assert sha(DATA/'training_arrays.npz')==lock['source_data_sha256']
    for name,expected in lock['code_sha256'].items():assert sha(WORK/name)==expected
    with np.load(DATA/'training_arrays.npz') as z:d={k:z[k] for k in z.files}
    with np.load(DATA/'measured_pairs.npz') as z:pairs={k:z[k] for k in z.files}
    train_mask=d['split']!=2
    if args.held_study=='none':test_mask=d['split']==2
    else:
        s=list(d['study_names']).index(args.held_study);groups=np.unique(d['group'][d['study']==s]);train_mask&=(d['study']!=s)&~np.isin(d['group'],groups);test_mask=d['study']==s
    tr,te=np.flatnonzero(train_mask),np.flatnonzero(test_mask);assert not set(d['group'][tr])&set(d['group'][te])
    med,support,weights,centers,quality,context=transformations(d,train_mask,True)
    nl=len(d['cell_line_names']);nc=len(med);parent=np.zeros(nl,dtype='int64');line_med=np.zeros(nl,dtype='float32');line_support=np.zeros(nl,dtype=bool)
    cell_sources=np.zeros(nc,dtype='int16');line_sources=np.zeros(nl,dtype='int16')
    for c in range(nc):cell_sources[c]=len(np.unique(d['study'][train_mask&(d['cell']==c)]))
    for line in range(nl):
        cells=np.unique(d['cell'][d['line']==line]);assert len(cells)<=1
        if len(cells):parent[line]=cells[0]
        ix=train_mask&(d['line']==line);line_support[line]=ix.any();line_med[line]=np.median(d['effect'][ix]) if ix.any() else med[parent[line]];line_sources[line]=len(np.unique(d['study'][ix]))
    h,local=global_bank(),local_bank();ab=avi_bank() if args.avi else None;transform=None
    if args.avi:
        values=ab[np.unique(d['bank_index'][tr])];known=values[:,2]==1;assert known.any();median=np.nanmedian(values[known,:2],axis=0);filled=np.where(np.isfinite(values[:,:2]),values[:,:2],median);mean,sd=filled.mean(0),filled.std(0);sd[sd<1e-5]=1
        transform={'median':median.tolist(),'mean':mean.tolist(),'sd':sd.tolist()}
    torch.manual_seed(args.seed);model=LineAwareEffectNetwork(med,line_med,parent,lock['recipe'],args.seed,args.avi)
    arm='AVI' if args.avi else 'sequence';name=f'fixed_{args.held_study}_{arm}_{args.seed}'
    directory=FINAL/'fits';directory.mkdir(exist_ok=True);record=directory/(name+'_record.json');assert not record.exists()
    checkpoint=directory/(name+'.safetensors');initial_hash=shared_initial_hash(model,args.avi)
    metadata={'recipe':lock['recipe'],'epochs':lock['epochs'],'seed':args.seed,'held_study':args.held_study,'arm':arm,
              'training_rows':len(tr),'prediction_rows':len(te),'cell_medians':med.tolist(),'line_medians':line_med.tolist(),
              'line_parent_cell':parent.tolist(),'supported_cells':support.tolist(),'supported_lines':line_support.tolist(),
              'cell_names':d['cell_names'].tolist(),'cell_line_names':d['cell_line_names'].tolist(),
              'training_cell_source_counts':cell_sources.tolist(),'training_line_source_counts':line_sources.tolist(),
              'activity_and_precision_transformations':context,'AVI_transform':transform,
              'shared_initial_parameter_sha256':initial_hash,'training_index_sha256':hashlib.sha256(tr.tobytes()).hexdigest(),
              'prediction_index_sha256':hashlib.sha256(te.tobytes()).hexdigest(),'settings_lock_sha256':sha(lock_path),
              'code_sha256':sha(Path(__file__)),'PID':os.getpid(),'fresh_model_performance_read':False}
    write_json(directory/(name+'_started.json'),metadata)
    optimizer=torch.optim.AdamW([p for p in model.parameters() if p.requires_grad],lr=3e-4,weight_decay=.001)
    tp={k:p[np.all(train_mask[p],axis=1)] for k,p in pairs.items()}
    def forward(ids):
        banks=d['bank_index'][ids];gg=torch.from_numpy(h[banks]);ll=torch.from_numpy(local[banks]);c=d['cell'][ids];line=d['line'][ids]
        cf=np.eye(nc,dtype='float32')[c];cf*=support[c,None];lf=np.eye(nl,dtype='float32')[line];lf*=line_support[line,None]
        av=None
        if args.avi:
            a=ab[banks];a[:,:2]=(np.where(np.isfinite(a[:,:2]),a[:,:2],median)-mean)/sd;av=torch.from_numpy(np.clip(a,-12,12))
        return model(gg[:,0],gg[:,1],ll[:,0],ll[:,1],torch.from_numpy(c.astype('int64')),torch.from_numpy(line.astype('int64')),torch.from_numpy(cf),torch.from_numpy(lf),av)
    rng=np.random.default_rng(args.seed);updates=0;curve=[];start=time.monotonic()
    for epoch in range(1,lock['epochs']+1):
        torch.manual_seed(args.seed+epoch);model.train();order=rng.permutation(tr);losses=[]
        for i in range(0,len(order),1024):
            main_ids=order[i:i+1024];chunks=[main_ids];counts=[]
            for key in ['local','cross']:
                p=tp[key];count=32 if len(p) else 0;counts.append(count)
                if count:chunks.append(p[rng.integers(len(p),size=count)].reshape(-1))
            ids=np.concatenate(chunks);optimizer.zero_grad();pred=forward(ids)
            batch={k:torch.from_numpy(d[k][ids]) for k in ['effect','signed_effect','activity_reference','activity_alternate']};batch['activity_center']=torch.from_numpy(centers[ids])
            loss=activity_loss(pred,batch,torch.from_numpy(weights[ids]),torch.from_numpy(quality[ids]),len(main_ids),*counts);assert torch.isfinite(loss);loss.backward();torch.nn.utils.clip_grad_norm_([p for p in model.parameters() if p.requires_grad],1.);optimizer.step();updates+=1;losses.append(float(loss.detach()))
        curve.append({'epoch':epoch,'mean_objective':float(np.mean(losses)),'elapsed_seconds':time.monotonic()-start});write_json(directory/(name+'_curve.json'),curve)
        print('FIXED V9',name,epoch,round(time.monotonic()-start,1),flush=True)
    model.eval();output={k:np.empty(len(te),dtype='float32') for k in ['magnitude','signed_effect','activity_reference','activity_alternate']}
    with torch.no_grad():
        for i in range(0,len(te),2048):
            got=forward(te[i:i+2048])
            for k in output:output[k][i:i+2048]=got[k].numpy()
    save_file({k:v.detach().contiguous() for k,v in model.state_dict().items()},checkpoint)
    predictions=directory/(name+'_predictions.npz');np.savez_compressed(predictions,observation_index=te,original_split=d['split'][te],**output)
    write_json(record,{**metadata,'actual_backprop_updates':updates,'checkpoint':str(checkpoint),'checkpoint_sha256':sha(checkpoint),
                       'predictions':str(predictions),'predictions_sha256':sha(predictions),'completed_epochs':lock['epochs'],
                       'no_calibration_or_test_checkpoint_selection':True,'scope':'Molecular prediction only; final untouched-study evaluation remains.'})
    print('FIXED V9 COMPLETE',name,flush=True)


if __name__=='__main__':main()

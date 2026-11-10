"""Real cell-matched reader summaries from complete source DNA; no new-test inputs."""
import sys,importlib.util,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
import torch
from safetensors.torch import load_file

PREP=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen/candidate_pretrained_readers'
DEST=ROOT/'outputs/sleep_frailty_ml_next'
CELLS=['K562','HepG2','WTC11']
NAMES=['reference_activity_readout','alternate_activity_readout','signed_readout_difference','absolute_readout_difference',
       'global_difference_mean_absolute','global_difference_RMS','global_cosine_distance',
       'variant_local_difference_mean_absolute','variant_local_difference_RMS','variant_local_cosine_distance',
       'global_context_mean','global_context_sd','variant_local_context_mean','variant_local_context_sd','global_difference_maximum_absolute','variant_local_difference_maximum_absolute']


def local_pool(a,positions):
    centers=torch.div(positions,16,rounding_mode='floor').clamp(0,a.shape[-1]-1);ids=centers[:,None]+torch.tensor([-1,0,1],device=a.device);valid=(ids>=0)&(ids<a.shape[-1]);ids=ids.clamp(0,a.shape[-1]-1)
    values=torch.gather(a,2,ids[:,None,:].expand(-1,a.shape[1],-1));return (values*valid[:,None,:]).sum(-1)/valid.sum(-1)[:,None]


def summarize(g,l,p,n):
    gr,ga=g[:n],g[n:];lr,la=l[:n],l[n:];pr,pa=p[:n],p[n:];gd=ga-gr;ld=la-lr
    def cosine(a,b):return 1-(a*b).sum(1)/np.maximum(np.sqrt((a*a).sum(1)*(b*b).sum(1)),1e-12)
    return np.column_stack([pr,pa,pa-pr,abs(pa-pr),abs(gd).mean(1),np.sqrt((gd*gd).mean(1)),cosine(gr,ga),abs(ld).mean(1),np.sqrt((ld*ld).mean(1)),cosine(lr,la),
                            (gr.mean(1)+ga.mean(1))*.5,(gr.std(1)+ga.std(1))*.5,(lr.mean(1)+la.mean(1))*.5,(lr.std(1)+la.std(1))*.5,abs(gd).max(1),abs(ld).max(1)]).astype('float32')


def main():
    assert_previous_sources();DEST.mkdir(exist_ok=True);receipt=PREP/'cell_reader_tensor_and_architecture_audit.json';checked=json.loads(receipt.read_text());assert checked['passed']
    with np.load(DATA/'training_arrays.npz') as z:used=np.unique(z['bank_index'])
    assert len(used)==255400;needed=set(map(int,used));pairs={}
    for source in [V7/'sequence_bank.csv.gz',DATA/'extra_sequence_bank.csv.gz']:
        for r in csv_rows(source):
            index=int(r['bank_id'])
            if index in needed:
                assert index not in pairs;assert len(r['reference'])==len(r['alternate']);changes=[i for i,(a,b) in enumerate(zip(r['reference'],r['alternate'])) if a!=b];assert changes==[int(r['changed_position'])]
                pairs[index]=(r['reference'],r['alternate'],changes[0])
    assert set(pairs)==needed
    feature_path=DEST/'cell_matched_native_development_features.npy';assert not feature_path.exists()
    features=np.lib.format.open_memmap(feature_path,mode='w+',dtype='float32',shape=(int(used.max())+1,3,16));features[:]=np.nan
    physical=np.full((int(used.max())+1,4),-1,dtype='int32');groups=collections.defaultdict(list)
    alphabet='AGCT';lut=np.full(256,255,dtype='uint8');lut[[ord(c) for c in alphabet]]=np.arange(4,dtype='uint8')
    for index,(a,b,j) in pairs.items():groups[len(a)].append(index);physical[index]=[len(a),j,alphabet.index(a[j]),alphabet.index(b[j])]
    device='mps' if torch.backends.mps.is_available() else 'cpu';torch.set_num_threads(2);arch=V8/'sources/legnet_model.py';assert sha(arch)==checked['author_architecture_sha256'];spec=importlib.util.spec_from_file_location('pinned_new_feature_architecture',arch);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    models={}
    for cell in CELLS:
        target=PREP/f'{cell}_author_test1_val2.safetensors';assert sha(target)==checked['checks'][cell]['tensor_file_sha256'];cfg=json.loads((PREP/f'{cell}_author_config.json').read_text());keys=['stem_ch','stem_ks','ef_ks','ef_block_sizes','pool_sizes','resize_factor'];model=module.LegNet(in_ch=4,**{k:cfg[k] for k in keys});model.load_state_dict(load_file(target));model.eval();models[cell]=model.to(device)
    plan={'purpose':'Additional real blood/liver/stem-cell pretrained features for development. No claim of stronger prediction.',
          'primary_source_data_sha256':sha(DATA/'training_arrays.npz'),'actual_physical_pairs':len(used),'reader_order':CELLS,'features_per_reader':NAMES,
          'sequence_policy':'Use the entire actual source reference/alternate DNA pair in its recorded assay orientation. Group by true length; no padding, invented flanks, length equalization, new outcome labels or source-activity input.',
          'length_counts':{str(length):len(ids) for length,ids in groups.items()},'physical_fields':['actual_DNA_length','actual_changed_index','reference_base_AGCT','alternate_base_AGCT'],
          'readers_frozen_during_feature_construction':True,'future_public_blood_or_other_library_inputs_used':False,'author_readout_unit':'Original reader activity feature, not a measured current-study activity or a disease-effect prediction.',
          'runtime':{'python':sys.version,'numpy':np.__version__,'torch':torch.__version__,'device':device},'PID':os.getpid()}
    write_json(DEST/'cell_matched_feature_plan.json',plan);start=time.monotonic();done=0;max_replay=0.;max_local=0.
    for length,ids in sorted(groups.items()):
        for offset in range(0,len(ids),256):
            batch=ids[offset:offset+256];rr=[pairs[i][0] for i in batch];aa=[pairs[i][1] for i in batch];positions=[pairs[i][2] for i in batch];tokens=np.stack([lut[np.frombuffer(s.encode(),dtype='uint8')] for s in rr+aa]);assert tokens.max()<4
            x=torch.nn.functional.one_hot(torch.from_numpy(tokens.astype('int64')),4).permute(0,2,1).to(device=device,dtype=torch.float32);pos=torch.tensor(positions+positions,device=device)
            for c,cell in enumerate(CELLS):
                model=models[cell]
                with torch.no_grad():
                    acts=model.main(model.stem(x));g=acts.mean(-1);local=local_pool(acts,pos);p=model.head(model.mapper(acts).mean(-1)).squeeze(-1).cpu().numpy()
                    if offset==0:
                        direct=model(x).cpu().numpy();assert direct.shape==p.shape;max_replay=max(max_replay,float(abs(direct-p).max()));assert max_replay<1e-4
                        v=acts.cpu().numpy();expected=[]
                        for row,j in zip(v,positions+positions):
                            center=min(j//16,row.shape[-1]-1);expected.append(row[:,max(0,center-1):min(row.shape[-1],center+2)].mean(1))
                        max_local=max(max_local,float(abs(np.array(expected)-local.cpu().numpy()).max()))
                values=summarize(g.cpu().numpy(),local.cpu().numpy(),p,len(batch));assert values.shape==(len(batch),16) and np.isfinite(values).all();features[batch,c]=values
            done+=len(batch)
            if offset%(256*32)==0:features.flush();print('CELL-MATCHED NATIVE DNA FEATURES',done,'/',len(used),'length',length,'seconds',round(time.monotonic()-start,1),flush=True)
    features.flush();assert np.isfinite(features[used]).all();assert max_replay<1e-4 and max_local<1e-4
    np.save(DEST/'native_development_physical_metadata.npy',physical)
    write_json(DEST/'cell_matched_feature_receipt.json',{'features_sha256':sha(feature_path),'physical_metadata_sha256':sha(DEST/'native_development_physical_metadata.npy'),'plan_sha256':sha(DEST/'cell_matched_feature_plan.json'),
              'physical_pairs':len(used),'all_used_features_finite':True,'independent_forward_replay_maximum_error':max_replay,'independent_NumPy_local_pool_maximum_error':max_local,
              'new_test_cases_or_effects_used':False,'new_model_trained':False,'elapsed_seconds':time.monotonic()-start})
    assert_previous_sources();print('CELL-MATCHED NATIVE DEVELOPMENT FEATURES COMPLETE',flush=True)


if __name__=='__main__':main()

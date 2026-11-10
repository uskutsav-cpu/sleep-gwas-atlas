"""Retain fixed reader-channel coordinates from complete original DNA pairs."""
import importlib.util,time
from channel_common import *
from extract_cell_matched_development_features import PREP,CELLS,local_pool,summarize
import torch
from safetensors.torch import load_file

def main():
    assert_frozen_previous();plan=get_plan();out=DEST/'projected_reader_channels.npy';assert not out.exists()
    with np.load(OLD/'coverage_training_arrays.npz') as z:used=np.unique(z['bank_index'])
    needed=set(map(int,used));pairs={}
    for source in [V7/'sequence_bank.csv.gz',DATA/'extra_sequence_bank.csv.gz',OLD/'OL27_development_sequence_bank.csv.gz']:
        for r in csv_rows(source):
            index=int(r['bank_id'])
            if index not in needed:continue
            assert index not in pairs;positions=[j for j,(a,b) in enumerate(zip(r['reference'],r['alternate'])) if a!=b];assert positions==[int(r['changed_position'])]
            assert len(r['reference'])==len(r['alternate']);pairs[index]=(r['reference'],r['alternate'],positions[0])
    assert pairs.keys()==needed
    n=int(used.max())+1;features=np.lib.format.open_memmap(out,mode='w+',dtype='float16',shape=(n,3,4,16));features[:]=np.nan
    projections=[]
    for c in range(3):
        rng=np.random.default_rng(plan['projection_seed']+c);q,_=np.linalg.qr(rng.standard_normal((128,16)));assert abs(q.T@q-np.eye(16)).max()<1e-12;projections.append(q.astype('float32'))
    projections=np.asarray(projections);np.save(DEST/'fixed_channel_projections.npy',projections)
    checked=json.loads((PREP/'cell_reader_tensor_and_architecture_audit.json').read_text());assert checked['passed'];arch=V8/'sources/legnet_model.py';assert sha(arch)==checked['author_architecture_sha256'];spec=importlib.util.spec_from_file_location('fixed_channel_author_architecture',arch);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    device='mps' if torch.backends.mps.is_available() else 'cpu';torch.set_num_threads(2);models={}
    for cell in CELLS:
        target=PREP/f'{cell}_author_test1_val2.safetensors';assert sha(target)==checked['checks'][cell]['tensor_file_sha256'];cfg=json.loads((PREP/f'{cell}_author_config.json').read_text());keys=['stem_ch','stem_ks','ef_ks','ef_block_sizes','pool_sizes','resize_factor'];model=module.LegNet(in_ch=4,**{k:cfg[k] for k in keys});model.load_state_dict(load_file(target));models[cell]=model.eval().to(device)
    oldbank=np.load(OLD/'cell_matched_native_development_features.npy',mmap_mode='r');extra=np.load(OLD/'coverage_appended_reader_features.npy',mmap_mode='r');offset=len(oldbank)
    alphabet='AGCT';lut=np.full(256,255,dtype='uint8');lut[[ord(c) for c in alphabet]]=np.arange(4,dtype='uint8');groups=collections.defaultdict(list)
    for index,(a,b,j) in pairs.items():groups[len(a)].append(index)
    start=time.monotonic();done=0;projection_error=0.;summary_error=0.;quantization_relative=0.;zeroed_nonzero=0;max_storage=0.
    for length,ids in sorted(groups.items()):
        for start_at in range(0,len(ids),256):
            batch=ids[start_at:start_at+256];rr=[pairs[i][0] for i in batch];aa=[pairs[i][1] for i in batch];positions=[pairs[i][2] for i in batch];tokens=np.stack([lut[np.frombuffer(s.encode(),dtype='uint8')] for s in rr+aa]);assert tokens.max()<4
            x=torch.nn.functional.one_hot(torch.from_numpy(tokens.astype('int64')),4).permute(0,2,1).to(device=device,dtype=torch.float32);pos=torch.tensor(positions+positions,device=device);nb=len(batch)
            for c,cell in enumerate(CELLS):
                model=models[cell]
                with torch.no_grad():
                    acts=model.main(model.stem(x));g=acts.mean(-1).cpu().numpy();local=local_pool(acts,pos).cpu().numpy();p=model.head(model.mapper(acts).mean(-1)).squeeze(-1).cpu().numpy()
                # Projection is fixed before labels and does not fit any new source.
                blocks=np.stack([(g[:nb]+g[nb:])*.5,g[nb:]-g[:nb],(local[:nb]+local[nb:])*.5,local[nb:]-local[:nb]],axis=1)
                projected=(blocks@projections[c]).astype('float32');projected[:,[1,3]]*=plan['difference_storage_multiplier'];stored=projected.astype('float16');assert np.isfinite(stored).all()
                features[batch,c]=stored;max_storage=max(max_storage,float(abs(projected).max()));scale=np.maximum(abs(projected),1e-3);quantization_relative=max(quantization_relative,float((abs(stored.astype('float32')-projected)/scale).max()));zeroed_nonzero+=int(np.sum((stored==0)&(projected!=0)))
                if start_at==0:
                    direct=np.empty_like(projected,dtype='float64')
                    for bi in range(nb):
                        for kind in range(4):direct[bi,kind]=np.array([sum(float(a)*float(b) for a,b in zip(blocks[bi,kind],projections[c,:,j])) for j in range(16)])
                    direct[:,[1,3]]*=plan['difference_storage_multiplier'];projection_error=max(projection_error,float(abs(direct-projected).max()));expected=summarize(g,local,p,nb)
                    previous=np.array([oldbank[i,c] if i<offset else extra[i-offset,c] for i in batch]);summary_error=max(summary_error,float(abs(expected-previous).max()))
            done+=nb
            if start_at%(256*32)==0:
                features.flush();write_json(DEST/'channel_feature_progress.json',{'physical_pairs_done':done,'total':len(used),'elapsed_seconds':time.monotonic()-start,'PID':os.getpid()});print('CHANNEL FEATURES',done,'/',len(used),'length',length,'seconds',round(time.monotonic()-start,1),flush=True)
    features.flush();assert done==len(used);assert projection_error<.01 and summary_error<1e-3 and quantization_relative<.001
    # Inspect in blocks to avoid allocating a full second feature bank.
    for i in range(0,len(used),4096):assert np.isfinite(features[used[i:i+4096]]).all()
    write_json(DEST/'channel_feature_receipt.json',{'passed':True,'physical_pairs':len(used),'shape':list(features.shape),'features_sha256':sha(out),'projections_sha256':sha(DEST/'fixed_channel_projections.npy'),
              'plan_sha256':sha(DEST/'finite_channel_plan.json'),'source_data_sha256':plan['source_data_sha256'],'reader_order':CELLS,'coordinate_order':['global_mean','signed_global_difference_x1000','local_mean','signed_local_difference_x1000'],
              'independent_float64_projection_max_error':projection_error,'recomputed_scalar_summary_max_error':summary_error,'maximum_relative_quantization_error_with_floor_1e_minus3':quantization_relative,'nonzero_coordinates_rounded_to_zero':zeroed_nonzero,'largest_stored_magnitude':max_storage,
              'actual_length_counts':{str(k):len(v) for k,v in groups.items()},'pretrained_weights_updated':False,'previous_blood_or_OL33_cases_added':False,'no_padding_or_invented_context':True,
              'device':device,'code_sha256':sha(Path(__file__)),'elapsed_seconds':time.monotonic()-start});assert_frozen_previous();print('CHANNEL FEATURES COMPLETE',flush=True)

if __name__=='__main__':main()

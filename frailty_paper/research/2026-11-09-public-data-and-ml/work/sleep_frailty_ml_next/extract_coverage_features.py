"""Small appended feature bank; never rewrite original cached development features."""
import sys,importlib.util,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
from extract_cell_matched_development_features import PREP,DEST,CELLS,NAMES,local_pool,summarize
import torch
from safetensors.torch import load_file

def main():
    assert_previous_sources();receipt=json.loads((DEST/'coverage_expansion_receipt.json').read_text());assert receipt['passed']
    source=DEST/'OL27_development_sequence_bank.csv.gz';assert sha(source)==receipt['appended_sequence_bank_sha256']
    rows=list(csv_rows(source));n=len(rows);assert n==receipt['appended_physical_pairs'];start_bank=receipt['appended_bank_start']
    out=DEST/'coverage_appended_reader_features.npy';assert not out.exists();features=np.full((n,3,16),np.nan,dtype='float32');physical=np.full((n,4),-1,dtype='int32')
    checked=json.loads((PREP/'cell_reader_tensor_and_architecture_audit.json').read_text());assert checked['passed'];arch=V8/'sources/legnet_model.py';assert sha(arch)==checked['author_architecture_sha256']
    spec=importlib.util.spec_from_file_location('pinned_coverage_architecture',arch);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    device='mps' if torch.backends.mps.is_available() else 'cpu';torch.set_num_threads(2);models={}
    for cell in CELLS:
        target=PREP/f'{cell}_author_test1_val2.safetensors';assert sha(target)==checked['checks'][cell]['tensor_file_sha256'];cfg=json.loads((PREP/f'{cell}_author_config.json').read_text());keys=['stem_ch','stem_ks','ef_ks','ef_block_sizes','pool_sizes','resize_factor'];model=module.LegNet(in_ch=4,**{k:cfg[k] for k in keys});model.load_state_dict(load_file(target));model.eval();models[cell]=model.to(device)
    alphabet='AGCT';lut=np.full(256,255,dtype='uint8');lut[[ord(c) for c in alphabet]]=np.arange(4,dtype='uint8');start=time.monotonic();replay=0.;local_error=0.
    for i in range(0,n,256):
        batch=rows[i:i+256];positions=[int(r['changed_position']) for r in batch];rr=[r['reference'] for r in batch];aa=[r['alternate'] for r in batch]
        assert all(len(a)==len(b)==200 for a,b in zip(rr,aa));assert [int(r['bank_id']) for r in batch]==list(range(start_bank+i,start_bank+i+len(batch)))
        for j,(a,b,p) in enumerate(zip(rr,aa,positions)):
            assert [k for k,(x,y) in enumerate(zip(a,b)) if x!=y]==[p];physical[i+j]=[200,p,alphabet.index(a[p]),alphabet.index(b[p])]
        tokens=np.stack([lut[np.frombuffer(s.encode(),dtype='uint8')] for s in rr+aa]);assert tokens.max()<4
        x=torch.nn.functional.one_hot(torch.from_numpy(tokens.astype('int64')),4).permute(0,2,1).to(device=device,dtype=torch.float32);pos=torch.tensor(positions+positions,device=device)
        for c,cell in enumerate(CELLS):
            model=models[cell]
            with torch.no_grad():
                acts=model.main(model.stem(x));g=acts.mean(-1);local=local_pool(acts,pos);p=model.head(model.mapper(acts).mean(-1)).squeeze(-1).cpu().numpy()
                if i==0:
                    replay=max(replay,float(abs(model(x).cpu().numpy()-p).max()));v=acts.cpu().numpy();expected=[]
                    for row,j in zip(v,positions+positions):
                        center=min(j//16,row.shape[-1]-1);expected.append(row[:,max(0,center-1):min(row.shape[-1],center+2)].mean(1))
                    local_error=max(local_error,float(abs(np.array(expected)-local.cpu().numpy()).max()))
            features[i:i+len(batch),c]=summarize(g.cpu().numpy(),local.cpu().numpy(),p,len(batch))
        if i%(256*16)==0:print('APPENDED DEVELOPMENT FEATURES',i+len(batch),'/',n,round(time.monotonic()-start,1),'seconds',flush=True)
    assert np.isfinite(features).all() and replay<1e-4 and local_error<1e-4
    np.save(out,features);np.save(DEST/'coverage_appended_physical_metadata.npy',physical)
    write_json(DEST/'coverage_feature_receipt.json',{'passed':True,'appended_pairs':n,'appended_bank_start':start_bank,'reader_order':CELLS,'features_per_reader':NAMES,
              'features_sha256':sha(out),'physical_metadata_sha256':sha(DEST/'coverage_appended_physical_metadata.npy'),'source_sequence_bank_sha256':sha(source),
              'independent_forward_replay_maximum_error':replay,'independent_local_pool_maximum_error':local_error,'device':device,
              'actual_full_DNA_length':200,'padded_or_invented_DNA':False,'new_test_inputs_or_outcomes_used':False,'pretrained_weights_updated':False,
              'code_sha256':sha(Path(__file__)),'elapsed_seconds':time.monotonic()-start});assert_previous_sources();print('APPENDED DEVELOPMENT FEATURES COMPLETE',flush=True)

if __name__=='__main__':main()

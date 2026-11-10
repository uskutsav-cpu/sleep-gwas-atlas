"""Physical identities and full-DNA features only, with no reserved outcomes."""
import sys,importlib.util,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
from audit_and_lock_library_measurements import write_gzip_csv
from extract_cell_matched_development_features import PREP,DEST,CELLS,NAMES,local_pool,summarize
import torch
from safetensors.torch import load_file
SCREEN=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen'

def main():
    assert_previous_sources();out=DEST/'reserved_reader_features.npy';assert not out.exists()
    assert json.loads((DEST/'library_measurement_role_lock.json').read_text())['passed']
    roles={(r['library'],r['original_variant_GRCh37']):r for r in csv_rows(SCREEN/'Siraj_unused_libraries/pretraining_guarded_source_roles.csv')}
    available={(r['library'],r['cell_file'],r['original_variant_GRCh37']):r for r in csv_rows(DEST/'library_measurement_availability_manifest.csv.gz')}
    physical=[];lookup={};observations=[]
    def add(design,metadata):
        key=(design['reference'],design['alternate'])
        if key not in lookup:
            lookup[key]=len(physical);physical.append({'bank_id':len(physical),'reference':key[0],'alternate':key[1],'changed_position':int(design['changed_position']),
                                                      'variant_GRCh38':design['variant_GRCh38'],'chromosome':design['chromosome'],'position_hg38':int(design['position_hg38'])})
        observations.append({'row_id':len(observations),'bank_index':lookup[key],**metadata})
    for design in csv_rows(SCREEN/'Siraj_unused_libraries/full_reference_and_shared_core_roles.csv'):
        role=roles[(design['library'],design['original_variant_GRCh37'])]
        if role['role']!='reserved_primary':continue
        for cell,line,group in [('HepG2','HepG2','hepatocyte_tumor'),('K562','K562','hematopoietic'),('SKNSH','SK-N-SH','neuronal_tumor')]:
            a=available[(design['library'],cell,design['original_variant_GRCh37'])]
            add(design,{'source':'Siraj2026','library':design['library'],'original_variant_GRCh37':design['original_variant_GRCh37'],
                        'source_SNP':'','comparison':'','window':'','source_cell':cell,'cell_group':group,'cell_line':line,'exact_line_known':True,
                        'component':'Siraj|'+role['component_root'],'variant_GRCh38':design['variant_GRCh38'],
                        'coefficient_available':a['finite_original_magnitude_available'],'primary_scope':'Held library from same publication/lab; OL27-connected reserved rows descriptive only.'})
    for design in csv_rows(SCREEN/'candidate_strict_no_published_pretraining_neighborhood.csv'):
        for cell in ['K562','HL60','CHRF','THP1']:
            add(design,{'source':'Tardaguila2026','library':'blood148','original_variant_GRCh37':design['original_GRCh37_variant'],
                        'source_SNP':design['source_SNP'],'comparison':design['comparison'],'window':design['window'],'source_cell':cell,
                        'cell_group':'hematopoietic','cell_line':cell,'exact_line_known':cell=='K562','component':'Tardaguila|'+design['component_root'],
                        'variant_GRCh38':design['variant_GRCh38'],'coefficient_available':True,
                        'primary_scope':'Independent publication; exact known K562 primary. Unknown blood lines diagnostic and withheld.'})
    write_gzip_csv(DEST/'reserved_physical_sequence_bank.csv.gz',physical);write_gzip_csv(DEST/'reserved_observation_metadata.csv.gz',observations)
    # Same full-DNA feature pipeline as development; no learned head or source labels.
    checked=json.loads((PREP/'cell_reader_tensor_and_architecture_audit.json').read_text());assert checked['passed'];arch=V8/'sources/legnet_model.py';assert sha(arch)==checked['author_architecture_sha256']
    spec=importlib.util.spec_from_file_location('pinned_reserved_architecture',arch);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    device='mps' if torch.backends.mps.is_available() else 'cpu';torch.set_num_threads(2);models={}
    for cell in CELLS:
        path=PREP/f'{cell}_author_test1_val2.safetensors';assert sha(path)==checked['checks'][cell]['tensor_file_sha256'];cfg=json.loads((PREP/f'{cell}_author_config.json').read_text());keys=['stem_ch','stem_ks','ef_ks','ef_block_sizes','pool_sizes','resize_factor'];model=module.LegNet(in_ch=4,**{k:cfg[k] for k in keys});model.load_state_dict(load_file(path));models[cell]=model.eval().to(device)
    features=np.full((len(physical),3,16),np.nan,dtype='float32');pp=np.full((len(physical),4),-1,dtype='int32');alphabet='AGCT';lut=np.full(256,255,dtype='uint8');lut[[ord(c) for c in alphabet]]=np.arange(4,dtype='uint8');by_length=collections.defaultdict(list)
    for i,r in enumerate(physical):
        a,b=r['reference'],r['alternate'];p=int(r['changed_position']);assert [j for j,(x,y) in enumerate(zip(a,b)) if x!=y]==[p]
        assert len(a)==len(b);pp[i]=[len(a),p,alphabet.index(a[p]),alphabet.index(b[p])];by_length[len(a)].append(i)
    start=time.monotonic();replay=0.;local_error=0.
    for length,ids in sorted(by_length.items()):
        for start_at in range(0,len(ids),256):
            batch=ids[start_at:start_at+256];rr=[physical[i]['reference'] for i in batch];aa=[physical[i]['alternate'] for i in batch];pos=[physical[i]['changed_position'] for i in batch]
            tokens=np.stack([lut[np.frombuffer(s.encode(),dtype='uint8')] for s in rr+aa]);assert tokens.max()<4;x=torch.nn.functional.one_hot(torch.from_numpy(tokens.astype('int64')),4).permute(0,2,1).to(device=device,dtype=torch.float32);positions=torch.tensor(pos+pos,device=device)
            for c,cell in enumerate(CELLS):
                model=models[cell]
                with torch.no_grad():
                    acts=model.main(model.stem(x));g=acts.mean(-1);local=local_pool(acts,positions);p=model.head(model.mapper(acts).mean(-1)).squeeze(-1).cpu().numpy()
                    if start_at==0:
                        replay=max(replay,float(abs(model(x).cpu().numpy()-p).max()));expected=[]
                        for row,j in zip(acts.cpu().numpy(),pos+pos):
                            center=min(j//16,row.shape[-1]-1);expected.append(row[:,max(0,center-1):min(row.shape[-1],center+2)].mean(1))
                        local_error=max(local_error,float(abs(np.array(expected)-local.cpu().numpy()).max()))
                features[batch,c]=summarize(g.cpu().numpy(),local.cpu().numpy(),p,len(batch))
            if start_at%(256*16)==0:print('RESERVED DNA FEATURES',length,start_at+len(batch),'/',len(ids),flush=True)
    assert np.isfinite(features).all() and replay<1e-4 and local_error<1e-4;np.save(out,features);np.save(DEST/'reserved_physical_metadata.npy',pp)
    write_json(DEST/'reserved_input_receipt.json',{'passed':True,'observations':len(observations),'physical_pairs':len(physical),'actual_lengths':{str(k):len(v) for k,v in by_length.items()},
              'observations_by_source_library':dict(collections.Counter(r['source']+'|'+r['library'] for r in observations)),
              'features_sha256':sha(out),'physical_metadata_sha256':sha(DEST/'reserved_physical_metadata.npy'),
              'observation_metadata_sha256':sha(DEST/'reserved_observation_metadata.csv.gz'),'sequence_bank_sha256':sha(DEST/'reserved_physical_sequence_bank.csv.gz'),
              'magnitude_model_fitted_or_run':False,'reserved_numeric_outcomes_joined':False,'unfamiliar_cell_lines_will_be_withheld':True,
              'independent_forward_replay_maximum_error':replay,'independent_local_pool_maximum_error':local_error,
              'code_sha256':sha(Path(__file__)),'elapsed_seconds':time.monotonic()-start});assert_previous_sources();print('RESERVED PHYSICAL INPUTS COMPLETE',flush=True)

if __name__=='__main__':main()

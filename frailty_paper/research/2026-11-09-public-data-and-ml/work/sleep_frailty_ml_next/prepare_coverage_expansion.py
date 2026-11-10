"""One locked OL27 expansion; old arrays and partitions remain exact prefixes."""
import sys,time,itertools
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
from audit_and_lock_library_measurements import partition,write_gzip_csv

DEST=ROOT/'outputs/sleep_frailty_ml_next'

def main():
    start=time.monotonic();assert_previous_sources();receipt=json.loads((DEST/'library_measurement_role_lock.json').read_text());assert receipt['passed']
    source=DEST/'OL27_development_measurements.csv.gz';assert sha(source)==receipt['development_measurements_sha256']
    output=DEST/'coverage_training_arrays.npz';assert not output.exists()
    with np.load(DATA/'training_arrays.npz') as z:old={k:z[k] for k in z.files}
    with np.load(DATA/'measured_pairs.npz') as z:old_pairs={k:z[k] for k in z.files}
    rows=list(csv_rows(source));rows.sort(key=lambda r:(r['original_variant_GRCh37'],r['line']))
    old_n=len(old['effect']);bank_start=len(np.load(DEST/'cell_matched_native_development_features.npy',mmap_mode='r'))
    names={k:list(old[k]) for k in old if k.endswith('_names')};indices={k:{v:i for i,v in enumerate(values)} for k,values in names.items()};new=[];physical=[];bank_lookup={};bank_group={};components={};position={}
    def code(key,value):
        if value not in indices[key]:indices[key][value]=len(names[key]);names[key].append(value)
        return indices[key][value]
    for r in rows:
        assert r['library']=='OL27' and r['role']=='development_candidate'
        component=r['component_root'];split=partition(component);assert split==int(r['development_split'])
        group=code('locus_names','OL27comp|'+component);components.setdefault(component,split);assert components[component]==split
        identity=(r['reference'],r['alternate'])
        if identity not in bank_lookup:
            bank_lookup[identity]=bank_start+len(physical);physical.append({k:r[k] for k in ['original_variant_GRCh37','variant_GRCh38','reference','alternate','changed_position','chromosome','position_hg38','component_root']})
        bank=bank_lookup[identity];bank_group.setdefault(bank,group);assert bank_group[bank]==group;position[bank]=int(r['position_hg38'])
        new.append({'bank_index':bank,'effect':float(r['effect']),'cell':code('cell_names',r['cell']),'line':code('cell_line_names',r['line']),
                    'study':code('study_names','Siraj'),'experiment':code('experiment_names','Siraj|OL27|'+r['line']),
                    'context':code('context_names','Siraj|'+r['cell']+'|absolute_log2_effect'),
                    'activity_reference':np.nan,'activity_alternate':np.nan,'signed_effect':np.nan,
                    'measurement_standard_error':float(r['measurement_standard_error']),
                    'standard_error_kind':2 if np.isfinite(float(r['measurement_standard_error'])) else 0,
                    'original_V7_observation_index':-1,'original_V8_observation_index':-1,'source_previously_external_in_V7':0,'group':group,'split':split})
    assert len(names['study_names'])==len(old['study_names'])
    assert len(names['cell_names'])==len(old['cell_names']) and len(names['cell_line_names'])==len(old['cell_line_names'])
    combined={}
    for key,value in old.items():
        if key.endswith('_names'):combined[key]=np.array(names[key]);continue
        extra=np.array([r[key] for r in new],dtype=value.dtype);combined[key]=np.concatenate([value,extra]);assert np.array_equal(combined[key][:old_n],value,equal_nan=True)
    new_ids=np.arange(old_n,len(combined['effect']),dtype='int32');cross=[];local=[]
    grouping=collections.defaultdict(list)
    for idx in new_ids:grouping[int(combined['bank_index'][idx])].append(int(idx))
    for ids in grouping.values():
        for a,b in itertools.combinations(sorted(ids,key=lambda j:(combined['cell'][j],combined['line'][j])),2):
            if combined['cell'][a]!=combined['cell'][b]:cross.append((a,b))
    order=sorted(map(int,new_ids),key=lambda j:(combined['experiment'][j],combined['group'][j],position[int(combined['bank_index'][j])]))
    for a,b in zip(order,order[1:]):
        ba,bb=map(int,combined['bank_index'][[a,b]]);gap=position[bb]-position[ba]
        if combined['experiment'][a]==combined['experiment'][b] and combined['group'][a]==combined['group'][b] and 0<=gap<=25 and ba!=bb:local.append((a,b))
    pairs={key:np.concatenate([old_pairs[key],np.array(value,dtype='int32').reshape(-1,2)]) for key,value in [('local',local),('cross',cross)]}
    for key,p in pairs.items():
        assert np.array_equal(p[:len(old_pairs[key])],old_pairs[key]);assert np.all(combined['split'][p[:,0]]==combined['split'][p[:,1]]) and np.all(combined['group'][p[:,0]]==combined['group'][p[:,1]])
    for i in range(3):
        for j in range(i+1,3):assert not set(combined['group'][combined['split']==i])&set(combined['group'][combined['split']==j])
    for i,r in enumerate(physical):r['bank_id']=bank_start+i
    np.savez_compressed(output,**combined);np.savez_compressed(DEST/'coverage_measured_pairs.npz',**pairs);write_gzip_csv(DEST/'OL27_development_sequence_bank.csv.gz',physical)
    result={'passed':True,'original_observations':old_n,'added_observations':len(new),'total_observations':len(combined['effect']),
            'appended_physical_pairs':len(physical),'appended_bank_start':bank_start,'all_original_array_prefixes_equal_including_NaNs':True,
            'all_original_group_and_partition_assignments_unchanged':True,'all_original_measured_pair_prefixes_equal':True,
            'added_pairs':{'local':len(local),'cross':len(cross)},'total_pairs':{k:len(v) for k,v in pairs.items()},
            'new_studies_added':0,'no_reserved_source_components_or_test_targets_added':True,'source_role_lock_sha256':sha(DEST/'library_measurement_role_lock.json'),
            'original_training_arrays_sha256':sha(DATA/'training_arrays.npz'),'original_pairs_sha256':sha(DATA/'measured_pairs.npz'),
            'training_arrays_sha256':sha(output),'pairs_sha256':sha(DEST/'coverage_measured_pairs.npz'),'appended_sequence_bank_sha256':sha(DEST/'OL27_development_sequence_bank.csv.gz'),
            'signed_or_allele_activity_targets_inferred':False,'code_sha256':sha(Path(__file__)),'elapsed_seconds':time.monotonic()-start}
    write_json(DEST/'coverage_expansion_receipt.json',result);assert_previous_sources();print(json.dumps(result),flush=True)

if __name__=='__main__':main()

"""Independent direct checks of the augmented data, roles and held-out DNA."""
import bisect
from v9_common import *


def main():
    report=json.loads((OUT/'revealed_source_augmentation_audit.json').read_text())
    for name,h in report['outputs'].items():assert sha(DATA/name)==h,name
    with np.load(DATA/'training_arrays.npz') as z:d={k:z[k] for k in z.files}
    with np.load(V8/'measured_training_arrays_with_activity_and_SE.npz') as z:old={k:z[k] for k in z.files}
    with np.load(DATA/'partition_bank_metadata.npz') as z:p={k:z[k] for k in z.files}
    b0=report['old_feature_bank_offset'];present=d['original_V8_observation_index']>=0;ii=d['original_V8_observation_index'][present]
    for key in ['effect','cell','line','study','experiment','context','activity_reference','activity_alternate','signed_effect','measurement_standard_error','standard_error_kind']:
        assert np.array_equal(d[key][present],old[key][ii],equal_nan=True),key
    assert np.all(d['split'][present][old['split'][ii]==2]==2)
    assert np.all(d['split'][present][old['split'][ii]!=2]!=2)
    roles={int(i):int(s) for i,s in zip(d['bank_index'],d['split'])};groups={int(i):int(g) for i,g in zip(d['bank_index'],d['group'])}
    roots={int(i):int(r) for i,r in zip(p['bank_index'],p['component_root'])}
    order=np.lexsort((p['position'],p['chromosome']));a,b=order[:-1],order[1:];near=(p['chromosome'][a]==p['chromosome'][b])&(p['position'][b]-p['position'][a]<=2000)
    assert np.all(p['component_root'][a[near]]==p['component_root'][b[near]])
    accepted={int(i) for i,excluded in zip(p['bank_index'],p['excluded']) if not excluded}
    assert accepted==set(groups)
    for root in np.unique(p['component_root']):
        ids=p['bank_index'][p['component_root']==root]
        flags=p['excluded'][p['component_root']==root];assert np.all(flags==flags[0])
        if not flags[0]:assert len({roles[int(i)] for i in ids})==1 and len({groups[int(i)] for i in ids})==1
    fresh_points=collections.defaultdict(list);fresh_cores=set()
    for r in csv_rows(V8/'reserved_matoba_physical_design.csv'):
        fresh_points[r['chromosome']].append(int(r['position_hg38']))
        for s in [r['reference'],r['alternate']]:fresh_cores.update([s,rc(s)])
    for c in fresh_points:fresh_points[c].sort()
    seen={};checked=0;core_windows=0
    def inspect(row):
        nonlocal checked,core_windows
        i=int(row['bank_id'])
        if i not in accepted:return
        c=row['chromosome'];pos=int(row['position_hg38']);points=fresh_points[c];j=bisect.bisect_left(points,pos)
        assert not any(0<=k<len(points) and abs(points[k]-pos)<=2000 for k in [j-1,j])
        a,b=row['reference'],row['alternate'];p0=int(row['changed_position']);start=max(0,min(p0-75,len(a)-150))
        for s in [a,b]:
            for k in range(max(0,p0-149),min(p0,len(s)-150)+1):
                assert s[k:k+150] not in fresh_cores;core_windows+=1
            core=s[start:start+150];canonical=min(core,rc(core))
            if canonical in seen:
                other=seen[canonical];assert groups[i]==groups[other] and roles[i]==roles[other]
            else:seen[canonical]=i
        checked+=1
    for r in csv_rows(V7/'sequence_bank.csv.gz'):inspect(r)
    for r in csv_rows(DATA/'extra_sequence_bank.csv.gz'):inspect(r)
    source_coverage={}
    for c in np.unique(d['cell']):
        source_coverage[str(d['cell_names'][c])]=[str(d['study_names'][s]) for s in np.unique(d['study'][d['cell']==c])]
    write_json(OUT/'independent_augmented_partition_audit.json',{'passed':True,'observations_checked':len(d['effect']),'physical_pairs_checked':checked,'variant_containing150bp_windows_checked':core_windows,'all_existing_targets_activities_SE_unchanged':True,'old_calibration_roles_preserved_for_every_retained_observation':True,'whole_component_exclusions_verified':True,'all2kb_adjacencies_and_shared_canonical150bp_contexts_within_one_role':True,'fresh_test_2kb_and_actual150bp_overlaps':0,'cell_group_study_coverage':source_coverage,'source_data_sha256':sha(DATA/'training_arrays.npz'),'scope':'Independent data-leakage and partition audit. Earlier external sources are explicitly development, not new replication or independent people.'})
    assert_previous_sources();print('AUGMENTED PARTITIONS INDEPENDENTLY VERIFIED',len(d['effect']),checked,flush=True)


if __name__=='__main__':main()

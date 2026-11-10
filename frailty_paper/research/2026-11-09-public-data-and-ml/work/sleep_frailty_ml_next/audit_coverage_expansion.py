"""Separate reconstruction of development labels, identities, and prefix invariants."""
import sys,time,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
DEST=ROOT/'outputs/sleep_frailty_ml_next'
SOURCE=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen/Siraj_unused_libraries'

def main():
    start=time.monotonic();assert_previous_sources();r=json.loads((DEST/'coverage_expansion_receipt.json').read_text());assert r['passed']
    with np.load(DATA/'training_arrays.npz') as z:old={k:z[k] for k in z.files}
    with np.load(DEST/'coverage_training_arrays.npz') as z:d={k:z[k] for k in z.files}
    n=len(old['effect'])
    for key,value in old.items():
        prefix=d[key][:len(value)]
        assert np.array_equal(prefix,value,equal_nan=True) if value.dtype.kind in 'fc' else np.array_equal(prefix,value)
    roles={};reserved=set()
    for row in csv_rows(SOURCE/'pretraining_guarded_source_roles.csv'):
        roles[(row['library'],row['original_variant_GRCh37'])]=row
        if row['role'].startswith('reserved'):reserved.add(row['component_root'])
    sequence_rows=list(csv_rows(DEST/'OL27_development_sequence_bank.csv.gz'));seq={int(row['bank_id']):row for row in sequence_rows};pairkeys={(row['reference'],row['alternate']) for row in sequence_rows};literal_hits=0
    for path in [V7/'sequence_bank.csv.gz',DATA/'extra_sequence_bank.csv.gz']:
        for row in csv_rows(path):
            a,b=row['reference'],row['alternate']
            if (a,b) in pairkeys or (b,a) in pairkeys or (rc(a),rc(b)) in pairkeys or (rc(b),rc(a)) in pairkeys:literal_hits+=1
    assert literal_hits==0
    coefficients={}
    retrieval=json.loads((SOURCE/'targeted_retrieval_receipt.json').read_text())
    for item in retrieval['files']:
        if item['library']!='OL27':continue
        assert sha(ROOT/item['path'])==item['sha256']
        with gzip.open(ROOT/item['path'],'rt') as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                if row['comb']==row['SNP']+'_center_fwd_ref':coefficients[(item['cell'],'chr'+row['SNP'].removeprefix('chr'))]=row
    line_to_file={'HepG2':'HepG2','K562':'K562','SK-N-SH':'SKNSH'};seen=set();split_by_component={};split_groups={i:set(d['group'][d['split']==i]) for i in range(3)}
    for i in range(3):
        for j in range(i+1,3):assert not split_groups[i]&split_groups[j]
    for idx in range(n,len(d['effect'])):
        s=seq[int(d['bank_index'][idx])];key=('OL27',s['original_variant_GRCh37']);role=roles[key];assert role['role']=='development_candidate' and role['component_root'] not in reserved
        line=str(d['cell_line_names'][d['line'][idx]]);raw=coefficients[(line_to_file[line],s['original_variant_GRCh37'])]
        assert d['effect'][idx]==np.float32(abs(float(raw['Log2Skew'])))
        se=float(raw['Skew_SE']);assert d['measurement_standard_error'][idx]==np.float32(se)
        assert np.isnan(d['signed_effect'][idx]) and np.isnan(d['activity_reference'][idx]) and np.isnan(d['activity_alternate'][idx])
        assert str(d['study_names'][d['study'][idx]])=='Siraj';assert str(d['experiment_names'][d['experiment'][idx]])=='Siraj|OL27|'+line
        label=str(d['locus_names'][d['group'][idx]]);assert label=='OL27comp|'+role['component_root']
        value=int(hashlib.sha256(('OL27-coverage-20261108|'+role['component_root']).encode()).hexdigest()[:16],16)%1000000
        expected=0 if value<600000 else 1 if value<800000 else 2;assert int(d['split'][idx])==expected
        split_by_component.setdefault(role['component_root'],expected);assert split_by_component[role['component_root']]==expected
        identity=(int(d['bank_index'][idx]),line);assert identity not in seen;seen.add(identity)
    with np.load(DEST/'coverage_measured_pairs.npz') as z:pairs={k:z[k] for k in z.files}
    with np.load(DATA/'measured_pairs.npz') as z:
        for key in z.files:assert np.array_equal(pairs[key][:len(z[key])],z[key])
    for key,p in pairs.items():
        assert np.all(d['split'][p[:,0]]==d['split'][p[:,1]]) and np.all(d['group'][p[:,0]]==d['group'][p[:,1]])
        new=p[np.all(p>=n,axis=1)]
        if key=='cross':assert np.all(d['bank_index'][new[:,0]]==d['bank_index'][new[:,1]]) and np.all(d['cell'][new[:,0]]!=d['cell'][new[:,1]])
        else:
            assert np.all(d['experiment'][new[:,0]]==d['experiment'][new[:,1]])
            for a,b in new:assert 0<=int(seq[int(d['bank_index'][b])]['position_hg38'])-int(seq[int(d['bank_index'][a])]['position_hg38'])<=25
    result={'passed':True,'independently_reconstructed_original_coefficients':len(d['effect'])-n,'old_array_and_pair_prefixes_unchanged':True,
            'new_source_identity_and_component_partition_verified':True,'new_training_banks_matching_old_literal_pair_either_orientation':literal_hits,
            'original_signed_or_activity_targets_not_assumed':True,'whole_Siraj_omission_will_include_old_and_new_libraries':True,
            'coverage_arrays_sha256':sha(DEST/'coverage_training_arrays.npz'),'pairs_sha256':sha(DEST/'coverage_measured_pairs.npz'),
            'role_lock_sha256':sha(DEST/'library_measurement_role_lock.json'),'code_sha256':sha(Path(__file__)),
            'scope':'Independent target, ID, exact-pair and partition reconstruction. Wider history/shared-core and pretraining screens have separate source-hashed audits. No reserved model outcome evaluated.',
            'elapsed_seconds':time.monotonic()-start}
    write_json(DEST/'coverage_independent_data_audit.json',result);assert_previous_sources();print(json.dumps(result),flush=True)

if __name__=='__main__':main()

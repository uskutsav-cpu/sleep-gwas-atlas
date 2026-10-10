#!/usr/bin/env python3
"""Independent metadata binding audit of sealed acquisition v2; no launches."""
import datetime
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
P=ROOT/'sleep_unified_research_v4'
SEEN={}


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):
            h.update(b)
    return h.hexdigest()


def read(path):
    SEEN[str(path)]={'sha256':sha(path),'bytes':path.stat().st_size}
    return path.read_text()


def main():
    read(Path(__file__))
    pp=P/'manifests/extension_raw_acquisition_plan_v4_2.json'
    plan=json.loads(read(pp));assert sha(pp)=='b13346e198dc0d2f04d324aa0d1b76af33507ce72a299d9e3beccca8639ad544'
    old=json.loads(read(P/'manifests/extension_raw_acquisition_plan_v4.json'))
    source=Path(plan['executor_path']);body=read(source)
    snapshot=P/'reviews/independent_48_original_source_snapshot_v4_2.txt';assert not snapshot.exists();snapshot.write_text(body)
    proofs=[]
    for a,b in zip(old['members'],plan['members']):
        assert all(b[k]==v for k,v in a.items() if k!='body_path')
        assert Path(b['body_path']).name==a['filename'] and Path(b['body_path']).parent.name=='raw'
        r=json.loads(read(Path(b['historical_streaming_receipt'])));v=r['source_verification']
        checks={'same_locked_member':True,'original_sha':b['expected_sha256']==v['observed_sha256'],
                'original_md5':b['expected_md5']==v['expected_md5']==v['observed_md5'],
                'original_exact_size':b['expected_bytes']==v['expected_size_bytes']==v['observed_size_bytes'],
                'original_version':b['expected_s3_version_id']==v['http_version_id'],
                'original_etag':b['expected_etag']==v['http_etag'].strip('"'),
                'original_versioned_source':b['url']==v['source']==r['source_url'],
                'historical_verification_pass':v['verification_status']=='PASS'}
        proofs.append({'trait_id':b['extension_trait_id'],'checks':checks,'all_pass':all(checks.values())})
    assert len(proofs)==100 and all(r['all_pass'] for r in proofs)
    for s,h in plan['bound_sources'].items():
        q=Path(s);SEEN[s]={'sha256':sha(q),'bytes':q.stat().st_size};assert SEEN[s]['sha256']==h
    assert sum(m['expected_bytes'] for m in plan['members'])==plan['compressed_network_bytes']==227610388647
    resume=plan['explicit_resume'];partial=Path(resume['original_partial_path'])
    assert partial.stat().st_size==resume['prefix_bytes']==1129316352
    assert resume['source_trait']==plan['members'][0]['extension_trait_id']
    assert plan['additional_preserved_prefix_bytes']==resume['prefix_bytes']
    controls_path=P/'logs/extension_acquisition_fault_controls_v4_2.json'
    controls=json.loads(read(controls_path));assert controls['plan_sha256']==sha(pp) and controls['executor_sha256']==sha(source)
    assert controls['all_pass'] is True and len(controls['checks'])==5 and all(r['pass'] for r in controls['checks'])
    assert controls['network_requests']==controls['native_estimators_launched']==0
    test=P/'scripts/49_verify_acquisition_fault_controls.py';assert sha(test)==controls['test_script_sha256'];read(test)
    unchanged=all(sha(Path(s))==m['sha256'] for s,m in SEEN.items());assert unchanged
    receipt={'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'V2_IDENTITY_BINDINGS_PASS_FINAL_GUARD_AND_QUARANTINE_AMENDMENT_REQUIRED',
             'plan_sha256':sha(pp),'source_sha256':sha(source),'source_snapshot_sha256':sha(snapshot),'locked100_original_byte_MD5_SHA_version_ETag_bindings_pass':True,
             'member_checks':proofs,'bound_dependency_count':len(plan['bound_sources']),'all_bound_dependency_hashes_pass':True,'fault_control_receipt_binding_verified':True,
             'fault_controls_count':5,'fault_controls_not_independently_reexecuted':True,'network_or_native_launches_by_reviewer':0,
             'first_preserved_prefix_bytes':resume['prefix_bytes'],'prefix_content_independently_rehashed':False,
             'prefix_scope':'Stat size and frozen plan/preserved prior-receipt bindings verified;1.129GB prefix not rehashed again by reviewer.',
             'required_amendments':[
                 'Final post-cleanup and post-final-partial-hash resource/size/deadline/identity gates for each source and family.',
                 'Unconditional quarantine before any lock release even if receipt persistence raises or is interrupted.',
                 'Guard quarantine diagnostic IO and route catchable termination through owned cleanup; default SIGTERM currently bypasses finally.',
                 'Compare second retained-partial SHA/bytes to original expectations before rename, not merely record the values.',
                 'Reject empty or incomplete independent-review bindings in root admission.'],
             'size_cap_resume_limit':'Current curl max-filesize uses full expected size for resumed remaining transfer; sampled full-partial size and final full SHA remain mandatory. Confirm remaining-byte hard-cap semantics when revising.',
             'previous_review_qualification':'Earlier original45 review did not inspect historical streaming receipts; its SHA-unavailable wording should be read as absent from original45 plan, not absent from repository. All100 original SHA values are now independently confirmed available in bound historical receipts.',
             'full_gzip_integrity_or_harmonization_claimed':False,'new_scientific_result':False,'GWAS_body_opened':False,'all_consumed_hashes_unchanged':unchanged,'inputs':SEEN,'hash_buffer_bytes':65536}
    out=P/'reviews/independent_extension_acquisition_v2_receipt_v4_2.json';assert not out.exists();out.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['inputs','member_checks']},indent=2))


if __name__=='__main__':main()

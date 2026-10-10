#!/usr/bin/env python3
"""Metadata/code binding audit of exact real preprocessing plan; no asset parse."""
import hashlib
import json
from pathlib import Path
import os
import time

P=Path(__file__).resolve().parents[1];R=P/'reviews'
PLAN=P/'manifests/finngen_preprocessing_plan_v4_3_7.json'
EXPECTED='d75fb3f96795f38080081be89da10be04c82d74aa67d11821a7eac0d6e7a8046'

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def regular(p):
    p=Path(p)
    assert p.is_file() and not p.is_symlink() and not any(x.is_symlink() for x in p.parents)
    return p

def main():
    started=time.monotonic();assert sha(regular(PLAN))==EXPECTED
    p=json.loads(PLAN.read_text());d=p['dependencies_sha256'];base=json.loads(Path(p['baseline_gate_plan']).read_text())
    assert sha(p['baseline_gate_plan'])==p['baseline_gate_plan_sha256']=='05c0781aea2c7fe8f897d891e67bb741122e778382f35c8913946cf7469676b4'
    reference=set([p['reference'],p['chain'],p['allele_list']])
    reference.update(k for k in d if '/ref/eur_w_ld_chr/' in k)
    checked={};skipped={}
    for path,digest in d.items():
        regular(path)
        if path in reference:
            skipped[path]=dict(expected_sha256=digest,bytes=Path(path).stat().st_size,
                binding_method='Exact pinned reference hash inherited from original baseline or frozen reference identity; no new reference-body read.')
            if path in base['dependencies_sha256']:assert digest==base['dependencies_sha256'][path]
            continue
        assert sha(path)==digest
        checked[path]=digest
    frozen_source=dict(bytes=809346932,md5='f16c21acbf8b9ebfeb0f26a19f0ecfe3',sha256='353129dc8252461a2a97087ce6de39548d0b0bf1e488a6ef36da75e9f4690049')
    assert p['source_identity']==frozen_source and p['input_sha256']=={p['source']:frozen_source['sha256']}
    raw=regular(p['source']);assert raw.stat().st_size==frozen_source['bytes']
    acquired_path=P/'logs/finngen_R13_insomnia_acquisition_verified_v4_2.json'
    acquired=json.loads(acquired_path.read_text())
    assert acquired['source_path']==str(raw) and acquired['actual_sha256']==frozen_source['sha256'] and acquired['actual_md5']==frozen_source['md5'] and acquired['actual_size']==frozen_source['bytes']
    assert p['source_generation']=='1777989563097164'
    expected_refs={p['chain']:'14a712e8e147d9fc8e9d87d51977b46f6f8ddb93efbe5d0843d86b6205f587b1',
        p['reference']:'e6e4814d99a1eff91875014fe70c61f760182efdc1711ba6e155963cf5aa11f8',
        p['allele_list']:'ec73fca0b696e8beba465b51e52911676fcade375bcc9475d99c0ec30509d3ed'}
    assert all(d[x]==h for x,h in expected_refs.items()) and p['chain_bytes']==1246411
    assert p['stage']=='preprocessing' and len(p['jobs'])==1 and p['new_rg_commands']==0
    assert p['status']=='PREPARED_ONLY_NO_WORKER_LAUNCHED' and p['scientific_source_admitted'] is False and p['independent_replication_established'] is False
    assert p['all_190_original_commands_and_numerical_reviews_required'] is True
    g=p['guard'];assert (g['worker_count'],g['BLAS_threads'],g['internal_floor_bytes'],g['SSD_floor_bytes'],g['observed_aggregate_worker_RSS_limit_bytes'],g['new_output_limit_bytes'],g['deadline_seconds'])==(1,1,3<<30,5<<30,2<<30,2<<30,7200)
    assert p['global_reservation_bytes']==300<<30 and p['assumed_effective_N']==185146.70377332723
    ns=Path(p['namespace']);assert ns==Path(p['source_namespace'])/'pipeline_replay_v3_7'
    assert p['shared_heavy_worker_lock']==str(Path(p['source_namespace']).parents[1]/'native_heavy_worker.lock')
    assert p['derivative']==str(ns/'derived/insomnia.sumstats.gz') and p['preprocessing_receipt']==str(ns/'derived/insomnia.preprocessing.json')
    job=p['jobs'][0]
    expected_command=[base['python'],'-B',str(P/'scripts/56_preprocess_finngen_insomnia_feasibility_v3.py'),
        '--plan',str(PLAN),'--plan-sha256','{PLAN_SHA256}','--inherited-heavy-lock-fd','{HEAVY_LOCK_FD}']
    assert job['command_template']==expected_command and job['job_id']=='finngen_insomnia_preprocessing'
    assert job['result_receipt']==p['preprocessing_receipt'] and job['output_prefix']==str(ns/'derived/insomnia')
    assert job['worker_receipt']==str(ns/'receipts_v4/preprocessing.execution_receipt_v3_7.json')
    expected_reviews=[str(R/n) for n in ['independent_finngen_preprocessing_science_review_v4.md','independent_finngen_row_science_receipt_v4.json',
        'independent_finngen_operational_preflight_v4_2.md','independent_finngen_operational_binding_receipt_v4_2.json',
        'independent_finngen_operational_prelaunch_v4_3_7.md','independent_finngen_operational_prelaunch_v4_3_7.json','independent_finngen_operational_prelaunch_seal_v4_3_7.json']]
    assert p['required_independent_review_paths']==expected_reviews
    for review in expected_reviews[:4]:regular(review)
    stop_path=P/'logs/finngen_preprocessing_performance_stop_v4_3_2.json';stop=json.loads(stop_path.read_text())
    assert len(stop['preserved_artifact_sha256'])==9 and stop['owned_cleanup_verified'] is True and stop['remaining_owned_group_members']==[]
    assert stop['shared_heavy_mutex_release_verified'] is True and stop['completed_source_CRC_or_retained_row_proof'] is False
    old_checks={}
    for path,digest in stop['preserved_artifact_sha256'].items():
        assert d[path]==digest and sha(regular(path))==digest;old_checks[path]=digest
    old_master=json.loads((P/'logs/finngen_preprocessing_stage_receipt_v4_3_2.json').read_text())
    assert old_master['status']=='FAILED_PRESERVED' and old_master['owned_cleanup_verified'] is True
    old_plan=json.loads((P/'manifests/finngen_preprocessing_plan_v4_3_2.json').read_text())
    assert Path(old_plan['terminal_pending']).exists() and not Path(old_plan['terminal_seal']).exists()
    old_worker=json.loads(Path(old_plan['jobs'][0]['worker_receipt']).read_text())
    assert old_worker['owned_cleanup_verified'] is True and old_worker['process_group_teardown']['remaining_group_members']==[]
    outputs=[p['derivative'],p['preprocessing_receipt'],job['worker_receipt'],p['stage_receipt'],p['terminal_pending'],p['terminal_seal'],p['admission']]
    absent={x:not(Path(x).exists() or Path(x).is_symlink()) for x in outputs};assert all(absent.values())
    source_bound_only=dict(identity=frozen_source,generation=p['source_generation'],actual_file_regular_and_exact_size=True,
        acquisition_receipt_sha256=sha(acquired_path),raw_hash_not_recomputed_by_this_metadata_review=True,
        prospective_preparer_and_worker_full_hash_before_after_and_CRC_gates=True)
    receipt=dict(schema='independent_actual_finngen_preprocessing_plan_binding_v4_3_7',plan_path=str(PLAN),plan_sha256=EXPECTED,
        all_checks_pass=True,dependencies_count=len(d),fresh_metadata_code_runtime_sha256=checked,
        reference_body_hashes_bound_without_fresh_read=skipped,source_binding=source_bound_only,
        exact_single_preprocessing_command=expected_command,zero_fits=True,guard=g,global_reservation_bytes=p['global_reservation_bytes'],
        preserved9_failed_attempt_hashes_current=old_checks,old_failed_PENDING_and_empty_owned_group_preserved=True,
        new_attempt_outputs_and_admission_absent=absent,required7_review_paths=expected_reviews,
        original190_completion_addendum_sha256=d[str(P/'logs/native190_root_independent_completion_addendum_v4.json')],
        plan_unchanged=sha(PLAN)==EXPECTED,elapsed_seconds=time.monotonic()-started,
        source_or_reference_body_parse=False,actual_worker_fit_network_mutex=False)
    target=R/'independent_finngen_prelaunch_binding_receipt_v4_3_7.json'
    with target.open('x') as f:json.dump(receipt,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(all_pass=True,fresh_metadata_code_runtime=len(checked),reference_assets_bound_without_read=len(skipped),old_failure_artifacts=9)))

if __name__=='__main__':main()

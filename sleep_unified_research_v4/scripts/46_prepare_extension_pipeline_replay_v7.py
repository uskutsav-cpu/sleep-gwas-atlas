#!/usr/bin/env python3
"""Prepare additive worker-consumer correction; source identities remain deferred.

No source/reference/runtime payload is reread and no worker is launched here.
Future root admission and the unchanged executor gates must hash every input.
"""
import datetime
import json
from pathlib import Path
import shutil
from extension_replay_common_v5 import sha,write_new,physical_mount,acquisition_operational_gate

P=Path(__file__).resolve().parents[1]
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OUT=SSD/'extension_pipeline_replay_v7'


def main():
    prior=SSD/'extension_pipeline_replay_v6/extension_pipeline_replay_plan_v6.json'
    prior_sha='d0fd223d468f26583bfd2bd477a51548646b57e91275d7733eb1ae814b94b798'
    if sha(prior)!=prior_sha:raise RuntimeError('SEALED_V6_REPLAY_PLAN_CHANGED_BEFORE_PARSE')
    plan=json.loads(prior.read_text())
    reject=P/'reviews/extension_pipeline_v6_worker_consumer_addendum_v1_seal.json'
    reject_sha='35254ae81812f01121255ab34d7ce967d24e402b8a1b973007454c81a28c0d7c'
    if sha(reject)!=reject_sha:raise RuntimeError('V6_TARGETED_REJECTION_SEAL_CHANGED_BEFORE_PARSE')
    r=json.loads(reject.read_text())
    if r.get('verdict')!='MEASURED_CONSUMER_OMISSIONS_REJECT_UNLAUNCHED_V6_PRESERVED':raise RuntimeError('EXACT_V6_TARGETED_REJECTION_REQUIRED')
    review_map={**r['file_sha256'],**r['bound_prior_review_and_code_sha256'],**r['exact_teardown_helper_sha256'],str(reject):reject_sha}
    for path,digest in review_map.items():
        if sha(path)!=digest:raise RuntimeError('PRESERVED_V6_REJECTION_OR_METADATA_CHANGED: '+path)
    old=prior.parent
    for member in plan['members']:
        for key,value in list(member.items()):
            if isinstance(value,str) and value.startswith(str(old)+'/'):
                member[key]=value.replace(str(old)+'/',str(OUT)+'/',1)
            elif key in ['harmonize_command','munge_command']:
                member[key]=[v.replace(str(old)+'/',str(OUT)+'/',1) if v.startswith(str(old)+'/') else v for v in value]
    inherited=dict(plan['dependencies_sha256'])
    plan['dependencies_sha256'].update(review_map)
    for code in [Path(__file__),P/'scripts/47_run_extension_pipeline_replay_v7.py']:
        plan['dependencies_sha256'][str(code)]=sha(code)
    plan['preserved_v6_plan']=str(prior);plan['preserved_v6_plan_sha256']=prior_sha
    plan['preserved_v6_targeted_rejection_seal']=str(reject);plan['preserved_v6_targeted_rejection_seal_sha256']=reject_sha
    plan['schema']='historical_extension_raw_pipeline_checkpoint_replay_plan_v7'
    plan['prepared_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    plan['preparation_identity_scope']='CURRENT_NEW_CODE_AND_ALL_TARGETED_V6_REJECTION/PRIOR_REVIEW_METADATA;INHERITED_PRIOR_DEPENDENCY_IDENTITIES_DEFERRED'
    plan['inherited_dependency_sha256_not_freshly_verified_at_preparation']=inherited
    plan['fresh_source_reference_runtime_body_verification_at_preparation']=False
    plan['actual_execution_admission_requires_full_current_dependencies_including_archived_inputs']=True
    plan['worker_receipt_consumer_contract']='EXACT400_KEYS_AND_COMMANDS;FROZEN_SHA_BEFORE_PARSE_AFTER_CONSUMPTION;STRICT_RETURN0_STOPNULL_PLANUNCHANGEDTRUE_METADATAERRORS_EMPTY_OWNEDTRUE_REAPED_CLEANUPERROR_EMPTY;FAILUREADDENDUM_VETO;CURRENT_REGULAR_OUTPUT_MAP'
    plan['master_persistence_contract']='INTENDED_SERIALIZED_MASTER_SHA_BEFORE_WRITE_AND_POSTWRITE_MATCH;TERMINAL_EXPECTS_FIXED_INTENT'
    plan['execution_preconditions']='Genuine complete full100 raw8 producer plus unchanged full190 completion; separately reviewed actual v7 plan and root admission; exact400 worker receipts and unchanged full decompressed-template comparison; all current inputs fully hash-verified at admission/execution; original guards; no automatic retry.'
    if len(plan['members'])!=100 or len({m['extension_trait_id'] for m in plan['members']})!=100 or plan['total_owned_commands']!=400 or plan['estimator_calls']!=0:raise RuntimeError('UNCHANGED100_400_ZERO_FIT_CONTRACT_REQUIRED')
    acquisition_operational_gate(plan)
    plan['physical_SSD_preflight']=physical_mount()
    plan['resource_preflight']={'internal_free_bytes':shutil.disk_usage('/System/Volumes/Data').free,'SSD_free_bytes':shutil.disk_usage(SSD).free}
    if plan['resource_preflight']['internal_free_bytes']<plan['guard']['internal_floor_bytes'] or plan['resource_preflight']['SSD_free_bytes']<plan['guard']['SSD_floor_bytes']:raise RuntimeError('UNCHANGED_RESOURCE_GUARD_FAILS')
    if sha(prior)!=prior_sha or sha(reject)!=reject_sha:raise RuntimeError('PRIOR_OR_REJECTION_CHANGED_AFTER_CONSUMPTION')
    for path,digest in review_map.items():
        if sha(path)!=digest:raise RuntimeError('REVIEW_METADATA_CHANGED_AFTER_CONSUMPTION')
    if OUT.exists():raise RuntimeError('PRIVATE_V7_NAMESPACE_ALREADY_EXISTS_NO_AUTOMATIC_RETRY')
    for folder in ['harmonized','munged','qc','receipts','logs_v4','receipts_v4','tmp','cache']:(OUT/folder).mkdir(parents=True,exist_ok=False)
    out=OUT/'extension_pipeline_replay_plan_v7.json';write_new(out,plan)
    print(json.dumps({'plan':str(out),'sha256':sha(out),'members':100,'commands':400,'estimator_calls':0,'workers_launched':0,'body_payloads_read':0,'inherited_dependencies_deferred':len(inherited)},indent=2))


if __name__=='__main__':main()

#!/usr/bin/env python3
"""Freeze mixed preserved/new raw checkpoints and durable terminal protocol.

This preparation does not launch workers or read raw source bodies. Historical
science commands change only their private output and exact raw-location paths.
"""
import datetime
import json
from pathlib import Path
import shutil
from extension_replay_common_v5 import sha,write_new,check_bindings,physical_mount,acquisition_operational_gate

P=Path(__file__).resolve().parents[1]
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OUT=SSD/'extension_pipeline_replay_v6'


def main():
    prior=SSD/'extension_pipeline_replay_v5/extension_pipeline_replay_plan_v5.json'
    prior_sha='766245a2fc3bef811fdd1217cb1fbff47175c7d518b4337900ca9e5fb17c7866'
    if sha(prior)!=prior_sha:raise RuntimeError('SEALED_V5_REPLAY_PLAN_CHANGED')
    plan=json.loads(prior.read_text());check_bindings(plan,include_archived=True)
    acquisition=P/'manifests/extension_raw_acquisition_plan_v4_8.json'
    acquisition_sha='77bf2da98b851ac1aa9d414b2b695fb04a65b0c4ca987f7d85db2097e9e7ef8b'
    if sha(acquisition)!=acquisition_sha:raise RuntimeError('REVIEWED_V8_ACQUISITION_PLAN_CHANGED')
    raw_plan=json.loads(acquisition.read_text())
    admission=P/'manifests/extension_raw_acquisition_admission_v4_8.json'
    root=json.loads(admission.read_text());executor=P/'scripts/52_acquire_extension_raw_sources_v8.py'
    if root.get('execution_admitted') is not True or root.get('plan_sha256')!=acquisition_sha or root.get('executor_sha256')!=sha(executor):
        raise RuntimeError('EXACT_V8_ACQUISITION_ROOT_ADMISSION_REQUIRED')
    prior_origin_identities=plan['acquisition_operational_identity_by_origin']
    old_folder=str(prior.parent);new_folder=str(OUT)
    for member,raw in zip(plan['members'],raw_plan['members']):
        old_raw=member['raw']
        if {k:v for k,v in member['acquisition_member'].items() if k!='body_path'}!={k:v for k,v in raw.items() if k!='body_path'}:
            raise RuntimeError('SCIENTIFIC_RAW_SOURCE_MEMBER_CHANGED')
        for key,value in list(member.items()):
            if isinstance(value,str) and value.startswith(old_folder+'/'):member[key]=value.replace(old_folder+'/',new_folder+'/',1)
            elif key in ['harmonize_command','munge_command']:
                member[key]=[v.replace(old_folder+'/',new_folder+'/',1) if v.startswith(old_folder+'/') else raw['body_path'] if v==old_raw else v for v in value]
        member['raw']=raw['body_path'];member['acquisition_member']=raw
        if str(raw['index']) in raw_plan['reused_checkpoints']:
            reuse=raw_plan['reused_checkpoints'][str(raw['index'])]
            member['acquisition_receipt']=reuse['receipt_path']
            member['acquisition_origin_plan_sha256']=reuse['original_acquisition_plan_sha256']
        else:
            member['acquisition_receipt']=str(Path(raw_plan['pending_path']).parent/'receipts'/(raw['extension_trait_id']+'.json'))
            member['acquisition_origin_plan_sha256']=acquisition_sha
    new_origin_identity={'source_folder':str(Path(raw_plan['pending_path']).parent),'curl':'/usr/bin/curl',
        'per_body_seconds_limit':raw_plan['per_body_seconds_limit'],
        'resume_offset_by_index':{'5':raw_plan['explicit_resume']['prefix_bytes']}}
    reused_origins={m['acquisition_origin_plan_sha256'] for m in plan['members'] if str(m['index']) in raw_plan['reused_checkpoints']}
    if len(raw_plan['reused_checkpoints'])!=4 or set(raw_plan['reused_checkpoints'])!={'1','2','3','4'}:
        raise RuntimeError('FROZEN4_REUSED_ORIGINS_REQUIRED')
    plan['acquisition_operational_identity_by_origin']={origin:prior_origin_identities[origin] for origin in reused_origins}
    plan['acquisition_operational_identity_by_origin'][acquisition_sha]=new_origin_identity
    plan['acquisition_plan']=str(acquisition);plan['acquisition_plan_sha256']=acquisition_sha
    files=[prior,acquisition,admission,executor,Path(__file__)]
    files.extend(P/'reviews'/n for n in ['genomicsem_extension_pipeline_v5_preflight.md','genomicsem_extension_pipeline_v5_preflight.json','genomicsem_extension_pipeline_v5_preflight_seal.json','independent_historical_munge_text_v1.sha256','independent_historical_munge_text_review_v1.md'])
    files.extend(P/'scripts'/n for n in ['47_run_extension_pipeline_replay_v6.py','extension_replay_common_v5.py','extension_replay_validate_v6.py','sensitivity_executor_v4_4.py','terminal_commit_common_v2.py','native_stage_completion_v4_3.py'])
    files.extend(P/'reviews'/n for n in ['genomicsem_native_stage_completion_review_v4_3.md','genomicsem_native_stage_completion_review_v4_3.json','genomicsem_native_stage_completion_review_v4_4.md','genomicsem_native_stage_completion_review_v4_4.json'])
    deps=plan['dependencies_sha256'];deps.update(raw_plan['bound_sources'])
    files.extend(Path(p) for p in root['independent_review_artifact_sha256'])
    files.extend(P/'reviews'/n for n in ['genomicsem_extension_checkpoint_review_v3.md','genomicsem_extension_checkpoint_review_v3.json','genomicsem_extension_checkpoint_review_seal_v3.json','independent_core_numerical_adjudication_v4.sha256','independent_whole_extension_v4.sha256','independent_whole_validation_v4.sha256'])
    for file in files:deps[str(file)]=sha(file)
    operation={str(admission):sha(admission),str(executor):sha(executor),str(acquisition):acquisition_sha}
    operation.update(root['independent_review_artifact_sha256'])
    plan['acquisition_execution_binding']={'root_admission':str(admission),'executor':str(executor),
        'executor_sha256':sha(executor),'independent_review_sha256':root['independent_review_artifact_sha256'],'sha256':operation}
    primary=P/'logs/extension_raw_acquisition_family_receipt_v4_8.json'
    plan['acquisition_family_receipt']=str(primary)
    plan['acquisition_family_terminal_contract']={'pending_path':raw_plan['pending_path'],
        'terminal_seal_path':raw_plan['terminal_seal_path'],
        'primary_receipt_paths':[str(primary),str(Path(raw_plan['pending_path']).parent/primary.name)],
        'reused_exact_source_count':4,'new_exact_source_count':96,
        'checkpoint_credit_while_pending':'INDIVIDUAL_EXACT_SUCCESSFUL_SOURCE_ONLY',
        'full_family_credit':'BOTH_HASH_BOUND_PRIMARY_COPIES_AND_FULL100_TERMINAL_SEAL_AND_ABSENT_PENDING_AND_FAILURE_ADDENDA'}
    plan['guard']['family_terminal_assembly_seconds']=600
    ledger=Path(raw_plan['global_resource_ledger_path'])
    if sha(ledger)!=raw_plan['global_resource_ledger_sha256']:raise RuntimeError('GLOBAL_SSD_RESERVATION_LEDGER_CHANGED')
    plan['global_resource_ledger_path']=str(ledger);plan['global_resource_ledger_sha256']=sha(ledger)
    plan['reservation_arithmetic']=json.loads(ledger.read_text())
    plan['preserved_v5_plan']=str(prior);plan['preserved_v5_plan_sha256']=prior_sha
    plan['schema']='historical_extension_raw_pipeline_checkpoint_replay_plan_v6'
    plan['prepared_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    plan['pipeline_terminal_protocol']='DURABLE_PENDING_BEFORE_WORK;EXACT_RESULT_RECEIPT_AND_TERMINAL_SEAL;POST_PERSISTENCE_IDENTITY_RESOURCE_TERMINATION_CHECKS;ABSENT_PENDING_AND_FAILURE_ADDENDA_REQUIRED'
    plan['execution_preconditions']='Reviewed full190 baseline and exact historical text-header compatibility seal; exact root admission; four original reused sources plus96 original v8 transfers; one shared heavy worker; full100 producer terminal and400 worker receipts/all consumed regular outputs; exact entire processed template; no automatic retry.'
    plan['physical_SSD_preflight']=physical_mount()
    plan['resource_preflight']={'internal_free_bytes':shutil.disk_usage('/System/Volumes/Data').free,'SSD_free_bytes':shutil.disk_usage(SSD).free}
    if plan['resource_preflight']['internal_free_bytes']<plan['guard']['internal_floor_bytes'] or plan['resource_preflight']['SSD_free_bytes']<plan['guard']['SSD_floor_bytes']:raise RuntimeError('UNCHANGED_RESOURCE_GUARD_FAILS')
    if len(plan['members'])!=100 or len({m['extension_trait_id'] for m in plan['members']})!=100:raise RuntimeError('FROZEN100_PANEL_CARDINALITY_DIFFERS')
    for folder in ['harmonized','munged','qc','receipts','logs_v4','receipts_v4','tmp','cache']:(OUT/folder).mkdir(parents=True,exist_ok=True)
    acquisition_operational_gate(plan);check_bindings(plan,include_archived=True)
    out=OUT/'extension_pipeline_replay_plan_v6.json';write_new(out,plan)
    print(json.dumps(dict(plan=str(out),sha256=sha(out),members=100,commands=plan['total_owned_commands'],workers_launched=0,raw_bodies_read=0),indent=2))


if __name__=='__main__':main()

#!/usr/bin/env python3
"""Freeze additive checkpoint execution ordering; no real workers or raw reads."""
import datetime
import json
import shutil
from pathlib import Path
from extension_replay_common_v2 import sha,write_new,check_bindings,physical_mount,acquisition_operational_gate

P=Path(__file__).resolve().parents[1]
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OUT=SSD/'extension_pipeline_replay_v2'


def main():
    prior=SSD/'extension_pipeline_replay_v1/extension_pipeline_replay_plan_v1.json'
    prior_sha='f90df481115b99ac0ea3bae9744f510faedba01f272dc1fb101f245d87fb2d91'
    if sha(prior)!=prior_sha:raise RuntimeError('SEALED_V1_REPLAY_PLAN_CHANGED')
    plan=json.loads(prior.read_text());check_bindings(plan,include_archived=True)
    old_folder=str(prior.parent);new_folder=str(OUT)
    for member in plan['members']:
        for key,value in list(member.items()):
            if isinstance(value,str) and value.startswith(old_folder+'/'):member[key]=value.replace(old_folder+'/',new_folder+'/',1)
            elif key in ['harmonize_command','munge_command']:member[key]=[v.replace(old_folder+'/',new_folder+'/',1) for v in value]
    deps=plan['dependencies_sha256']
    files=[prior,prior.parent/'preparation_controls_v1.json',prior.parent/'PREPARATION_HANDOFF_v1.json',prior.parent/'PREPARATION_HANDOFF_v1.sha256',OUT/'METHODS_REVIEW_v2.md',Path(__file__)]
    names=['47_run_extension_pipeline_replay_v2.py','extension_replay_common_v2.py','extension_replay_validate_v2.py','extension_replay_fault_controls_v2.py']
    files.extend(P/'scripts'/n for n in names)
    admission=P/'manifests/extension_raw_acquisition_admission_v4_4.json';original=json.loads(admission.read_text())
    executor=P/'scripts/52_acquire_extension_raw_sources_v4.py'
    if original['execution_admitted'] is not True or original['plan_sha256']!=plan['acquisition_plan_sha256'] or original['executor_sha256']!=sha(executor):raise RuntimeError('ACQUISITION_ROOT_ADMISSION_NOT_EXACT')
    files.extend([admission,executor]);files.extend(Path(v) for v in original['independent_review_artifact_sha256'])
    for file in files:deps[str(file)]=sha(file)
    operation={str(admission):sha(admission),str(executor):sha(executor),plan['acquisition_plan']:plan['acquisition_plan_sha256']}
    operation.update(original['independent_review_artifact_sha256'])
    plan['acquisition_execution_binding']={'root_admission':str(admission),'executor':str(executor),'executor_sha256':sha(executor),'independent_review_sha256':original['independent_review_artifact_sha256'],'sha256':operation}
    plan['preserved_v1_plan_sha256']=prior_sha;plan['preserved_v1_plan']=str(prior)
    plan['schema']='historical_extension_raw_pipeline_checkpoint_replay_plan_v2'
    plan['prepared_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    plan['checkpoint_receipt_binding_policy']='FIRST_IMMUTABLE_SUCCESSFUL_EXACT52_RECEIPT_SHA_BOUND_BEFORE_FIRST_SOURCE_COMMAND_AND_RECHECKED_AFTER_EACH_COMMAND_AND_FINAL_SEAL'
    plan['acquisition_family_receipt']=str(P/'logs/extension_raw_acquisition_family_receipt_v4_4.json')
    plan['guard']['checkpoint_poll_seconds']=5;plan['guard']['receipt_assembly_seconds']=30
    plan['execution_preconditions']='Root exact operational admission and reviewed full190 baseline before queue; each frozen source receives its own immutable successful52receipt/SHA binding before commands; unlocked bounded waiting for unavailable sources/heavy mutex; no automatic retry.'
    plan['initial_root_receipt_hash_map']='Optional subset of already completed immutable successful receipts; future receipts are SHA-bound individually at admitted checkpoint discovery, not outcome selection.'
    plan['physical_SSD_preflight']=physical_mount()
    plan['resource_preflight']={'internal_free_bytes':shutil.disk_usage('/System/Volumes/Data').free,'SSD_free_bytes':shutil.disk_usage(SSD).free}
    if plan['resource_preflight']['internal_free_bytes']<plan['guard']['internal_floor_bytes'] or plan['resource_preflight']['SSD_free_bytes']<plan['guard']['SSD_floor_bytes']:raise RuntimeError('ORIGINAL_RESOURCE_GUARD_FAILS')
    for folder in ['harmonized','munged','qc','receipts','logs_v4','receipts_v4','tmp','cache']:(OUT/folder).mkdir(parents=True,exist_ok=True)
    acquisition_operational_gate(plan);check_bindings(plan,include_archived=True)
    out=OUT/'extension_pipeline_replay_plan_v2.json';write_new(out,plan)
    print(json.dumps(dict(plan=str(out),sha256=sha(out),members=len(plan['members']),commands=plan['total_owned_commands'],workers_launched=0),indent=2))


if __name__=='__main__':main()

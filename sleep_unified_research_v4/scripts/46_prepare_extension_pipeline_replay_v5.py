#!/usr/bin/env python3
"""Prepare original extension replay with previously existing text header fix."""
import datetime
import json
from pathlib import Path
import shutil
from extension_replay_common_v4 import sha,write_new,check_bindings,physical_mount,acquisition_operational_gate

P=Path(__file__).resolve().parents[1]
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OUT=SSD/'extension_pipeline_replay_v5'


def main():
    prior=SSD/'extension_pipeline_replay_v4/extension_pipeline_replay_plan_v3.json'
    prior_sha='7c941a230287e04b42526b215b1f13373b850885fd9c8abacd2f0ab83731c690'
    if sha(prior)!=prior_sha:raise RuntimeError('PRESERVED_V4_EXTENSION_PLAN_CHANGED')
    plan=json.loads(prior.read_text());check_bindings(plan,include_archived=True)
    compatibility=P/'source_provenance/historical_munge_header_text_compatibility_v1.json'
    if sha(compatibility)!='8e9388c92d6801c299a1d61b7d59fa91da789fd2c2422f4f4983e721dbcdb41a':raise RuntimeError('HISTORICAL_MUNGE_COMPATIBILITY_CHANGED')
    compat=json.loads(compatibility.read_text());old_folder=str(prior.parent);new_folder=str(OUT)
    for member in plan['members']:
        for key,value in list(member.items()):
            if isinstance(value,str) and value.startswith(old_folder+'/'):member[key]=value.replace(old_folder+'/',new_folder+'/',1)
            elif key in ['harmonize_command','munge_command']:
                member[key]=[v.replace(old_folder+'/',new_folder+'/',1) if v.startswith(old_folder+'/') else v for v in value]
        original=compat['pinned_unmodified_munger']
        if member['munge_command'].count(original)!=1:raise RuntimeError('EXACT_ORIGINAL_STOCK_MUNGE_PATH_REQUIRED')
        member['munge_command']=[compat['candidate_executable'] if v==original else v for v in member['munge_command']]
    files=[prior,compatibility,Path(__file__),Path(compat['historical_metadata_path'])]
    files.extend(Path(path) for path in compat['source_and_copy_identity'])
    files.extend(P/'scripts'/n for n in ['47_run_extension_pipeline_replay_v5.py','extension_replay_validate_v5.py'])
    files.extend(P/'reviews'/n for n in ['genomicsem_extension_pipeline_v4_preflight.md','genomicsem_extension_pipeline_v4_preflight.json','genomicsem_extension_pipeline_v4_preflight_seal.json'])
    for file in files:plan['dependencies_sha256'][str(file)]=sha(file)
    plan.update(schema='historical_extension_raw_pipeline_checkpoint_replay_plan_v5',prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),preserved_v4_plan=str(prior),preserved_v4_plan_sha256=prior_sha,historical_munge_compatibility_provenance=str(compatibility),historical_munge_compatibility_provenance_sha256=sha(compatibility),historical_munge_qualification=compat['qualification'],execution_preconditions='Reviewed full190 baseline; exact independent historical-header compatibility qualifier and root admission; frozen3+97 v7 source checkpoints; one shared heavy worker; full100 terminal contract; full original serialized content match; no automatic retry.')
    if len(plan['members'])!=100 or plan['total_owned_commands']!=400:raise RuntimeError('EXACT_ORIGINAL100_COMMAND_PANEL_REQUIRED')
    plan['physical_SSD_preflight']=physical_mount()
    plan['resource_preflight']={'internal_free_bytes':shutil.disk_usage('/System/Volumes/Data').free,'SSD_free_bytes':shutil.disk_usage(SSD).free}
    if plan['resource_preflight']['internal_free_bytes']<plan['guard']['internal_floor_bytes'] or plan['resource_preflight']['SSD_free_bytes']<plan['guard']['SSD_floor_bytes']:raise RuntimeError('UNCHANGED_EXTENSION_RESOURCE_FLOOR_FAILS')
    for folder in ['harmonized','munged','qc','receipts','logs_v4','receipts_v4','tmp','cache']:(OUT/folder).mkdir(parents=True,exist_ok=True)
    acquisition_operational_gate(plan);check_bindings(plan,include_archived=True)
    out=OUT/'extension_pipeline_replay_plan_v5.json';write_new(out,plan)
    print(json.dumps({'plan':str(out),'sha256':sha(out),'members':100,'commands':400,'workers_launched':0,'raw_bodies_read':0}))


if __name__=='__main__':main()

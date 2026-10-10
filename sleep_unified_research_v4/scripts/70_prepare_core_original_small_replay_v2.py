#!/usr/bin/env python3
"""Prepare the original 10 small-input core harmonizations; no workers/raw reads."""
import datetime
import json
from pathlib import Path
import shutil
from extension_replay_common_v3 import sha,write_new,physical_mount,check_bindings

P=Path(__file__).resolve().parents[1]
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OUT=SSD/'core_pipeline/original_small_input_replay_v2'


def main():
    manifest=P/'reviews/independent_core_pipeline_replay_manifest_v4.json'
    if sha(manifest)!='451ea6d16263a467b8ac587e2069d66dfc7f736c9e694e23801efb5ba307413b':raise RuntimeError('SEALED45_CORE_DESIGN_CHANGED')
    design=json.loads(manifest.read_text())
    runtime=P/'source_provenance/core_declared_harmonization_runtime_candidate_v4.json'
    if sha(runtime)!='26014a05a9ad60d42e250e6a61ad47cfb3d01eacad0db84c3d888cc88e919b57':raise RuntimeError('QUALIFIED_RUNTIME_PROOF_CHANGED')
    rt=json.loads(runtime.read_text())
    baseline=SSD/'extension_pipeline_replay_v3/extension_pipeline_replay_plan_v3.json'
    if sha(baseline)!='565b34e997b2cab841f3a101cb4f8649c4c4660e56f5b1421ce4f68f29afe4ac':raise RuntimeError('BASELINE_BLUEPRINT_CHANGED')
    b=json.loads(baseline.read_text())
    compatibility=P/'source_provenance/historical_munge_header_text_compatibility_v1.json'
    if sha(compatibility)!='8e9388c92d6801c299a1d61b7d59fa91da789fd2c2422f4f4983e721dbcdb41a':raise RuntimeError('HISTORICAL_MUNGE_COMPATIBILITY_PROVENANCE_CHANGED')
    compat=json.loads(compatibility.read_text())
    prior_failure=SSD/'core_pipeline/original_small_input_replay_v1/core_original_small_execution_receipt_v1.json'
    failed=json.loads(prior_failure.read_text())
    if failed['status']!='FAILED_PRESERVED_NO_AUTOMATIC_RETRY' or not failed['owned_cleanup_verified']:raise RuntimeError('PRIOR_CORE_FAILURE_NOT_PRESERVED_AND_CLEANED')
    keys=['baseline_support_package','baseline_execution_plan_sha256','baseline_dependency_sha256',
          'baseline_relocated_dependency_sha256','baseline_input_sha256','original_190_jobs',
          'core_precision_adjudication','python','ldsc_dir','reference_prefix','environment']
    plan={k:b[k] for k in keys};members=[]
    for r in design['rows']:
        if r['prefilter'] is None and r['trait_id']!='bmi':continue
        old=r['new_private_output_namespace'];new=str(OUT/r['trait_id'])
        commands={}
        for label,command in r['command_templates'].items():
            commands[label]=None if command is None else [rt['command'][0] if value=='{FROZEN_HARMONIZATION_PYTHON}' else value.replace(old+'/',new+'/',1) if value.startswith(old+'/') else value for value in command]
        original_munge=str(P.parent.parent/'ldsc-code/munge_sumstats.py')
        if commands['munge_stock'].count(original_munge)!=1:raise RuntimeError('EXACT_ORIGINAL_MUNGE_EXECUTABLE_REQUIRED')
        commands['munge_stock']=[compat['candidate_executable'] if value==original_munge else value for value in commands['munge_stock']]
        trait=r['trait_id'];folder=OUT/trait
        members.append({'trait_id':trait,'original_design':r,'prefilter_command':commands['prefilter'],
            'harmonize_command':commands['harmonize_original'],'munge_command':commands['munge_stock'],
            'harmonized':str(folder/'harmonized'/(trait+'.harmonized.tsv.gz')),
            'harmonization_qc':str(folder/'harmonized'/(trait+'.qc.txt')),
            'munged':str(folder/'munged'/(trait+'.sumstats.gz')),'munged_prefix':str(folder/'munged'/trait),
            'source_gate_receipt':str(folder/'receipts/source_gate.json'),
            'comparison_receipt':str(folder/'receipts/full_content_QC_comparison.json')})
    if len(members)!=10 or sum(m['prefilter_command'] is not None for m in members)!=6:raise RuntimeError('ORIGINAL_SMALL_INPUT_MEMBERSHIP_DIFFERS')
    dependencies={path:meta['sha256'] for path,meta in design['input_metadata_code_reference_hashes'].items()}
    files=[manifest,runtime,baseline,Path(__file__),P/'manifests/global_SSD_resource_reservation_v4_3.json']
    files.extend(P/'scripts'/name for name in ['71_run_core_original_small_replay_v2.py','72_validate_core_original_small_replay_v2.py','sensitivity_executor_v4_4.py','terminal_commit_common_v2.py','extension_replay_common_v3.py','native_stage_completion_v4_3.py','30_prepare_and_run_ssd_native_campaign.py'])
    files.extend(P/'reviews'/name for name in ['core_bounded_parser_methods_review_v1.md','core_bounded_parser_methods_review_v1.json','terminal_commit_common_review_v2.md','terminal_commit_common_review_v2.json','terminal_commit_common_review_seal_v2.json'])
    files.extend([compatibility,prior_failure,SSD/'core_pipeline/original_small_input_replay_v1/core_pending_v1.json'])
    files.extend(Path(path) for path in compat['source_and_copy_identity'])
    files.append(Path(compat['historical_metadata_path']))
    failed_worker=SSD/'core_pipeline/original_small_input_replay_v1/receipts_v4/ms__munge.worker.json'
    files.extend([failed_worker,SSD/'core_pipeline/original_small_input_replay_v1/logs_v4/ms__munge.worker.stdout.log'])
    for file in files:dependencies[str(file)]=sha(file)
    plan.update(schema='original10_small_input_core_raw_replay_plan_v2',prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        members=members,member_count=10,original_core_family_size=45,remaining_large_input_traits=35,
        scientific_membership_or_threshold_changes=False,automatic_retry=False,estimator_calls=0,
        archived_input_sha256={},dependencies_sha256=dependencies,
        qualified_runtime_receipt=str(runtime),qualified_runtime_receipt_sha256=sha(runtime),
        runtime_qualification=rt['qualification'],harmonization_python=rt['command'][0],
        historical_munge_compatibility_provenance=str(compatibility),historical_munge_compatibility_provenance_sha256=sha(compatibility),
        historical_munge_qualification=compat['qualification'],preserved_failed_v1_receipt=str(prior_failure),preserved_failed_v1_receipt_sha256=sha(prior_failure),
        source_body_reads_in_preparation=0,workers_launched=0,
        global_resource_ledger_path=str(files[4]),global_resource_ledger_sha256=sha(files[4]),
        guard={'worker_count':1,'BLAS_threads':1,'internal_floor_bytes':3*(1<<30),'SSD_floor_bytes':5*(1<<30),
            'observed_aggregate_worker_RSS_limit_bytes':2*(1<<30),'new_output_limit_bytes':16*(1<<30),
            'deadline_seconds':96*3600,'per_worker_deadline_seconds':2*3600,'poll_seconds':2,
            'checkpoint_poll_seconds':5,'global_reservation_bytes':300*(1<<30),
            'shared_heavy_worker_lock':str(SSD/'native_heavy_worker.lock')},
        execution_preconditions='Exact root admission and independent review; all190 completed and collated; unchanged original01/21 and previously existing historical header-text-compatible stockmunge; exact before/after source identities; runtime candidate freshly verified read-only; full ordered recovered harmonized/munged content and all original QC steps; no historical binary identity claim.',
        terminal_contract='DURABLE_PRIVATE_OWNED_PENDING;EXACT_RECEIPT_AND_SEAL;POST_PERSISTENCE_IDENTITY_RESOURCE_TERMINATION_CHECKS;ABSENT_PENDING_AND_FAILURE_ADDENDA_REQUIRED')
    check_bindings(plan)
    for member in members:
        for folder in ['prefilter','harmonized','munged','receipts']:(OUT/member['trait_id']/folder).mkdir(parents=True,exist_ok=True)
    for folder in ['receipts_v4','logs_v4','tmp','cache']:(OUT/folder).mkdir(parents=True,exist_ok=True)
    plan['physical_SSD_preflight']=physical_mount()
    plan['resource_preflight']={'internal_free_bytes':shutil.disk_usage('/System/Volumes/Data').free,'SSD_free_bytes':shutil.disk_usage(SSD).free}
    if plan['resource_preflight']['internal_free_bytes']<plan['guard']['internal_floor_bytes'] or plan['resource_preflight']['SSD_free_bytes']<plan['guard']['SSD_floor_bytes']:raise RuntimeError('UNCHANGED_RESOURCE_GUARD_FAILS')
    path=OUT/'core_original_small_replay_plan_v2.json';write_new(path,plan)
    print(json.dumps({'plan':str(path),'sha256':sha(path),'members':10,'prefilter_commands':6,'harmonization_commands':10,'stock_munge_commands':10,'workers_launched':0,'raw_source_body_reads':0},indent=2))


if __name__=='__main__':main()

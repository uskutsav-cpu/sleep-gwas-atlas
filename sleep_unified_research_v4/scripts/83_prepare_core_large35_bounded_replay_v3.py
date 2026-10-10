#!/usr/bin/env python3
"""Prepare 35 original large raw chains through independently checked bounded science."""
import datetime
import json
from pathlib import Path
import shutil
from extension_replay_common_v3 import sha,write_new,physical_mount,check_bindings

P=Path(__file__).resolve().parents[1]
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OUT=SSD/'core_pipeline/large35_bounded_replay_v3'


def main():
    manifest=P/'reviews/independent_core_pipeline_replay_manifest_v4.json'
    if sha(manifest)!='451ea6d16263a467b8ac587e2069d66dfc7f736c9e694e23801efb5ba307413b':raise RuntimeError('SEALED45_CORE_DESIGN_CHANGED')
    design=json.loads(manifest.read_text())
    component=P/'reviews/independent_column_reader_v4_seal_v1.json'
    adapter=P/'reviews/core_bounded_harmonizer_source_correction_review_seal_v3.json'
    if sha(component)!='c364c8102e2abe2ac318f5407d383d7fba6ce78f21c78c0fa943d5da5378b240':raise RuntimeError('REVIEWED_COLUMN4_COMPONENT_SEAL_CHANGED')
    if sha(adapter)!='1b2ce879da41c09cd0494cfa6b40fb7d16176c45c6b93ce918b3d84042798e50':raise RuntimeError('REVIEWED_BOUNDED3_ADAPTER_SEAL_CHANGED')
    component_review=json.loads(component.read_text());adapter_review=json.loads(adapter.read_text())
    if component_review['candidate_sha256']!='1e5a7f86ab51cd9cec2d705144ca7623bd8e57e218c3889cf1a94b4db7c5b375' or adapter_review['candidate_sha256']!='3b87e849e3001a2a2d30fc0971dd3bd97847736abb9bed0dff521d1b5e1d0210' or adapter_review['component4_sha256']!=component_review['candidate_sha256']:raise RuntimeError('EXACT_REVIEWED_COMPONENT_AND_ADAPTER_CANDIDATES_REQUIRED')
    for path,expected in [(P/'scripts/core_column_reader_v4.py',component_review['candidate_sha256']),(P/'scripts/core_bounded_harmonizer_v3.py',adapter_review['candidate_sha256'])]:
        if path.is_symlink() or sha(path)!=expected:raise RuntimeError('REVIEWED_CORE_CANDIDATE_CHANGED')
    for review in [component_review,adapter_review]:
        for path,artifact in review['artifacts'].items():
            if Path(path).is_symlink() or not Path(path).is_file() or sha(path)!=artifact['sha256']:raise RuntimeError('SEALED_COMPONENT_REVIEW_ARTIFACT_CHANGED')
    runtime=P/'source_provenance/core_declared_harmonization_runtime_candidate_v4.json'
    if sha(runtime)!='26014a05a9ad60d42e250e6a61ad47cfb3d01eacad0db84c3d888cc88e919b57':raise RuntimeError('QUALIFIED_RUNTIME_PROOF_CHANGED')
    rt=json.loads(runtime.read_text())
    baseline=SSD/'extension_pipeline_replay_v3/extension_pipeline_replay_plan_v3.json'
    if sha(baseline)!='565b34e997b2cab841f3a101cb4f8649c4c4660e56f5b1421ce4f68f29afe4ac':raise RuntimeError('BASELINE_BLUEPRINT_CHANGED')
    b=json.loads(baseline.read_text())
    compatibility=P/'source_provenance/historical_munge_header_text_compatibility_v1.json'
    if sha(compatibility)!='8e9388c92d6801c299a1d61b7d59fa91da789fd2c2422f4f4983e721dbcdb41a':raise RuntimeError('HISTORICAL_MUNGE_COMPATIBILITY_PROVENANCE_CHANGED')
    compat=json.loads(compatibility.read_text())
    keys=['baseline_support_package','baseline_execution_plan_sha256','baseline_dependency_sha256',
          'baseline_relocated_dependency_sha256','baseline_input_sha256','original_190_jobs',
          'core_precision_adjudication','python','ldsc_dir','reference_prefix','environment']
    plan={k:b[k] for k in keys};members=[]
    for r in design['rows']:
        if r['prefilter'] is not None or r['trait_id']=='bmi':continue
        old=r['new_private_output_namespace'];new=str(OUT/r['trait_id'])
        commands={}
        for label,command in r['command_templates'].items():
            commands[label]=None if command is None else [rt['command'][0] if value=='{FROZEN_HARMONIZATION_PYTHON}' else value.replace(old+'/',new+'/',1) if value.startswith(old+'/') else value for value in command]
        original_munge=str(P.parent.parent/'ldsc-code/munge_sumstats.py')
        if commands['munge_stock'].count(original_munge)!=1:raise RuntimeError('EXACT_ORIGINAL_MUNGE_EXECUTABLE_REQUIRED')
        commands['munge_stock']=[compat['candidate_executable'] if value==original_munge else value for value in commands['munge_stock']]
        trait=r['trait_id'];folder=OUT/trait
        original_command=commands['harmonize_original']
        if original_command[:3]!=[rt['command'][0],'-B',str(P.parent/'scripts/01_harmonize.py')]:raise RuntimeError('EXACT_ORIGINAL_HARMONIZER_COMMAND_REQUIRED')
        bounded=[rt['command'][0],'-B',str(P/'scripts/core_bounded_harmonizer_v3.py'),'--original',original_command[2],
            '--spool',str(folder/'ephemeral_bounded_spool'),'--expected-rows',str(r['harmonizer_input_rows']),
            '--expected-source-sha256',r['raw']['sealed_verified_sha256'],'--',*original_command[3:]]
        members.append({'trait_id':trait,'original_design':r,'prefilter_command':commands['prefilter'],
            'harmonize_command':bounded,'literal_original_harmonize_command':original_command,'ephemeral_spool':str(folder/'ephemeral_bounded_spool'),'bounded_receipt':str(folder/'receipts/bounded_harmonization.json'),'column_manifest':str(folder/'receipts/native_column_preparation.json'),'spool_inventory':str(folder/'receipts/ephemeral_spool_inventory.json'),'spool_cleanup_receipt':str(folder/'receipts/ephemeral_spool_cleanup.json'),'munge_command':commands['munge_stock'],
            'harmonized':str(folder/'harmonized'/(trait+'.harmonized.tsv.gz')),
            'harmonization_qc':str(folder/'harmonized'/(trait+'.qc.txt')),
            'munged':str(folder/'munged'/(trait+'.sumstats.gz')),'munged_prefix':str(folder/'munged'/trait),
            'source_gate_receipt':str(folder/'receipts/source_gate.json'),
            'comparison_receipt':str(folder/'receipts/full_content_QC_comparison.json')})
    if len(members)!=35 or any(m['prefilter_command'] is not None for m in members):raise RuntimeError('ORIGINAL_LARGE35_MEMBERSHIP_DIFFERS')
    dependencies={path:meta['sha256'] for path,meta in design['input_metadata_code_reference_hashes'].items()}
    files=[manifest,runtime,baseline,Path(__file__),P/'manifests/global_SSD_resource_reservation_v4_4.json']
    files.extend(P/'scripts'/name for name in ['84_run_core_large35_bounded_replay_v3.py','85_validate_core_large35_bounded_replay_v3.py','sensitivity_executor_v4_4.py','terminal_commit_common_v2.py','extension_replay_common_v3.py','native_stage_completion_v4_3.py','30_prepare_and_run_ssd_native_campaign.py'])
    files.extend(P/'reviews'/name for name in ['core_bounded_parser_methods_review_v1.md','core_bounded_parser_methods_review_v1.json','terminal_commit_common_review_v2.md','terminal_commit_common_review_v2.json','terminal_commit_common_review_seal_v2.json'])
    files.append(compatibility)
    files.extend(Path(path) for path in compat['source_and_copy_identity'])
    files.append(Path(compat['historical_metadata_path']))
    files.extend(P/'scripts'/name for name in ['core_bounded_harmonizer_v3.py','core_column_reader_v4.py','86_cleanup_core_bounded_spool_v3.py'])
    files.extend([component,adapter,P/'logs/core_small_v2_root_independent_adjudication_v1.json'])
    for seal_name,expected in [('independent_spool_disposal_seal_v1.json','642fa6e1fbd12d8e36bae5d7172db184ec3fed0851cfae9c31a6103a8e8e7217'),('core_large35_bounded_controller_review_seal_v1.json','06ba91909321a1e7d0f481bf565f5f0ad938dabd139995f803d1fe7ebdf098e1'),('independent_spool_disposal_seal_v2.json','0f892f5b5ed3897c689162b5f3a7f91ace5bf609cdfe578c59ef76c55f346fa7'),('core_large35_bounded_controller_review_seal_v2.json','9ef7aa0a983033cf1ae0cac6bc232f5804f6a16e64f3de3ae32cf56d0d15e12e')]:
        seal=P/'reviews'/seal_name
        if sha(seal)!=expected:raise RuntimeError('PRESERVED_PREDECESSOR_REJECTION_REVIEW_CHANGED')
        prior=json.loads(seal.read_text());files.append(seal)
        for path,artifact in prior.get('artifacts',prior.get('file_sha256',{})).items():
            digest=artifact['sha256'] if isinstance(artifact,dict) else artifact
            if Path(path).is_symlink() or sha(path)!=digest:raise RuntimeError('PREDECESSOR_REVIEW_ARTIFACT_CHANGED')
            files.append(Path(path))
    files.extend(Path(path) for review in [component_review,adapter_review] for path in review['artifacts'])
    for file in files:dependencies[str(file)]=sha(file)
    plan.update(schema='original35_large_input_core_bounded_raw_replay_plan_v3',prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        members=members,member_count=35,private_namespace=str(OUT),original_core_family_size=45,already_reproduced_small_input_traits=10,
        scientific_membership_or_threshold_changes=False,automatic_retry=False,estimator_calls=0,
        archived_input_sha256={},dependencies_sha256=dependencies,
        qualified_runtime_receipt=str(runtime),qualified_runtime_receipt_sha256=sha(runtime),
        runtime_qualification=rt['qualification'],harmonization_python=rt['command'][0],
        historical_munge_compatibility_provenance=str(compatibility),historical_munge_compatibility_provenance_sha256=sha(compatibility),
        historical_munge_qualification=compat['qualification'],
        ephemeral_spool_policy='DELETE_ONLY_OWN_GENERATED_SCRATCH_AFTER_EXACT_FULL_ORIGINAL_CONTENT_QC_PASS;FREEZE_COPIED_CANDIDATE_MANIFESTS_AND_ENTIRE_INVENTORY_BEFORE_DISPOSAL;PRESERVE_FAILED_SCRATCH_AND_ALL_ORIGINAL_INPUTS',
        source_body_reads_in_preparation=0,workers_launched=0,
        component_review_seal=str(component),component_review_seal_sha256=sha(component),
        bounded_adapter_review_seal=str(adapter),bounded_adapter_review_seal_sha256=sha(adapter),
        scientific_adapter_qualification='Exactly frozen35 original argv routes only; full native one-column/global sample conversion memory is not previously measured on real sources. Supervised source-specific attempts and exact full original content/QC pass remain mandatory. Components and tiny fixtures certify software behavior only.',
        global_resource_ledger_path=str(files[4]),global_resource_ledger_sha256=sha(files[4]),
        guard={'worker_count':1,'BLAS_threads':1,'internal_floor_bytes':3*(1<<30),'SSD_floor_bytes':5*(1<<30),
            'observed_aggregate_worker_RSS_limit_bytes':2*(1<<30),'new_output_limit_bytes':16*(1<<30),
            'deadline_seconds':96*3600,'per_worker_deadline_seconds':2*3600,'poll_seconds':2,
            'checkpoint_poll_seconds':5,'global_reservation_bytes':300*(1<<30),
            'shared_heavy_worker_lock':str(SSD/'native_heavy_worker.lock')},
        execution_preconditions='Exact root admission and independent review; all190 completed and collated; unchanged original01 scientific stages via separately reviewed bounded adapter and previously existing historical header-text-compatible stockmunge; exact before/after source identities; runtime candidate freshly verified read-only; full ordered recovered harmonized/munged content and all original QC steps; no historical binary identity claim.',
        terminal_contract='DURABLE_PRIVATE_OWNED_PENDING;EXACT_RECEIPT_AND_SEAL;POST_PERSISTENCE_IDENTITY_RESOURCE_TERMINATION_CHECKS;ABSENT_PENDING_AND_FAILURE_ADDENDA_REQUIRED')
    check_bindings(plan)
    for member in members:
        for folder in ['prefilter','harmonized','munged','receipts']:(OUT/member['trait_id']/folder).mkdir(parents=True,exist_ok=True)
    for folder in ['receipts_v4','logs_v4','tmp','cache']:(OUT/folder).mkdir(parents=True,exist_ok=True)
    plan['physical_SSD_preflight']=physical_mount()
    plan['resource_preflight']={'internal_free_bytes':shutil.disk_usage('/System/Volumes/Data').free,'SSD_free_bytes':shutil.disk_usage(SSD).free}
    if plan['resource_preflight']['internal_free_bytes']<plan['guard']['internal_floor_bytes'] or plan['resource_preflight']['SSD_free_bytes']<plan['guard']['SSD_floor_bytes']:raise RuntimeError('UNCHANGED_RESOURCE_GUARD_FAILS')
    path=OUT/'core_large35_bounded_replay_plan_v3.json';write_new(path,plan)
    print(json.dumps({'plan':str(path),'sha256':sha(path),'members':35,'prefilter_commands':0,'harmonization_commands':35,'stock_munge_commands':35,'source_checks':35,'full_comparisons':35,'owned_scratch_cleanup_commands':35,'workers_launched':0,'raw_source_body_reads':0},indent=2))


if __name__=='__main__':main()

"""Prepare only the actual empty-apnea-frame scientific correction continuation."""
import datetime,json
from pathlib import Path
from extension_replay_common_v3 import sha,write_new
from core_large35_checkpoint_adoption_v2 import regular,adoption_gate
P=Path(__file__).resolve().parents[1]
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OLD=SSD/'core_pipeline/large35_checkpoint_continuation_v1'
OUT=SSD/'core_pipeline/large35_checkpoint_continuation_v2'
OLD_SHA='716b9df426dc5291296a4dbbedd9f057f615168c8c61945da4ccb7dd074f7bb0'

def main():
    metadata={}
    def bind(path,expected=None):
        path=regular(path);digest=sha(path)
        if expected is not None and digest!=expected:raise RuntimeError('EXACT_SCIENTIFIC_CHECKPOINT_METADATA_CHANGED: '+str(path))
        metadata[str(path)]=digest;return digest
    old_path=OLD/'core_large35_checkpoint_continuation_plan_v1.json';bind(old_path,OLD_SHA)
    old=json.loads(old_path.read_text())
    failed_path=OLD/'core_large35_continuation_execution_receipt_v1.json';bind(failed_path);failed=json.loads(failed_path.read_text())
    if (failed['status']!='FAILED_PRESERVED_NO_AUTOMATIC_RETRY' or failed['plan_sha256']!=OLD_SHA
        or failed['owned_cleanup_verified'] is not True or failed['termination_requests']
        or failed['completed_members'] or len(failed['worker_receipt_sha256'])!=1
        or old['member_count']!=27 or old['members'][0]['trait_id']!='sleep_apnea'):
        raise RuntimeError('EXACT_ACTUAL_FIRST_APNEA_HARMONIZER_FAILURE_REQUIRED')
    bind(OLD/'core_large35_continuation_pending_v1.json')
    admission=P/'manifests/core_large35_checkpoint_continuation_admission_v1.json'
    bind(admission,failed['admission_sha256']);admit=json.loads(admission.read_text())
    if admit['execution_admitted'] is not True or admit['plan_sha256']!=OLD_SHA:
        raise RuntimeError('EXACT_OLD_ROOT_ADMISSION_AND_INDEPENDENT_REVIEW_REQUIRED')
    for path,digest in admit['independent_review_sha256'].items():bind(path,digest)
    bind(P/'reviews/core_large35_checkpoint_continuation_author_receipt_v1.json')
    member=old['members'][0];source_worker=OLD/'receipts_v4/sleep_apnea__source.worker.json'
    source_sha=bind(source_worker,failed['worker_receipt_sha256'][str(source_worker)])
    source_record=json.loads(source_worker.read_text())
    fail_worker=OLD/'receipts_v4/sleep_apnea__harmonize.worker.json';bind(fail_worker);record=json.loads(fail_worker.read_text())
    for r,status,rc in [(source_record,'WORKER_COMPLETE_VERIFIED',0),(record,'WORKER_FAILED_PRESERVED',1)]:
        if (r['status']!=status or r['returncode']!=rc or r['plan_sha256']!=OLD_SHA or r['stop_reason'] is not None
            or r['owned_cleanup_verified'] is not True or r['metadata_errors'] or r['process_group_teardown']['remaining_group_members']
            or r['process_group_teardown'].get('cleanup_error')):
            raise RuntimeError('EXACT_REAPED_SOURCE_AND_HARMONIZER_PRODUCERS_REQUIRED')
    if source_record['command']!=failed['expected_worker_commands'][str(source_worker)] or record['command']!=member['harmonize_command']:
        raise RuntimeError('EXACT_CURRENT_PRODUCER_COMMANDS_REQUIRED')
    for path,digest in {**source_record['output_sha256'],**record['output_sha256']}.items():bind(path,digest)
    bind(member['source_gate_receipt'],source_record['output_sha256'][member['source_gate_receipt']])
    source=json.loads(Path(member['source_gate_receipt']).read_text());raw=member['original_design']['raw']
    if source['status']!='EXACT_CORE_RAW_SOURCE_GATE_PASS' or source['trait']!='sleep_apnea' or source['plan_sha256']!=OLD_SHA or source['raw_sha256']!=raw['sealed_verified_sha256'] or source['raw_bytes']!=raw['bytes']:
        raise RuntimeError('EXACT_COMPLETED_CURRENT_APNEA_SOURCE_PROOF_REQUIRED')
    manifest=Path(member['ephemeral_spool'])/'columns/column_preparation_manifest.json';manifest_sha=bind(manifest)
    columns=json.loads(manifest.read_text())
    if (columns['source']!=raw['resolved_path'] or columns['source_sha256_before']!=raw['sealed_verified_sha256']
        or columns['source_sha256_after']!=raw['sealed_verified_sha256'] or columns['expected_source_sha256']!=raw['sealed_verified_sha256']
        or columns['expected_rows']!=20170208 or columns['expected_rows']!=member['original_design']['harmonizer_input_rows']
        or set(columns['typed_pieces'])!={'SNP','CHR','BP','A1','A2','FRQ','BETA','SE','P'}
        or any(len(items)!=404 for items in columns['typed_pieces'].values())
        or columns['coordinate_label_parsed'] or any(columns['preinference_token_pieces'].values())
        or columns['physical_shape_gate']['raw_gzip_EOF_and_CRC_reached'] is not True):
        raise RuntimeError('EXACT_COMPLETED_APNEA_COLUMN_CHECKPOINT_REQUIRED')
    checkpoint={'manifest_path':str(manifest),'manifest_sha256':manifest_sha,'producer_plan':str(old_path),'producer_plan_sha256':OLD_SHA,
        'source_gate_receipt':member['source_gate_receipt'],'source_gate_receipt_sha256':metadata[member['source_gate_receipt']],
        'failed_harmonizer_worker':str(fail_worker),'failed_harmonizer_worker_sha256':metadata[str(fail_worker)],
        'raw_sha256':raw['sealed_verified_sha256'],'expected_rows':20170208,'typed_column_count':9,'pieces_per_column':404,
        'scope':'ONLY_HASH_BOUND_COMPLETED_NATIVE_TYPED_COLUMNS;OLD_POST_LIFTOVER_AND_ORDINARY_CACHES_UNSEALED_NOT_ADMITTED',
        'piece_bytes_freshly_verified_in_preparation':False,'old_failed_spool_preserved':True}
    def relocate(value):
        if isinstance(value,dict):return {k:relocate(v) for k,v in value.items()}
        if isinstance(value,list):return [relocate(v) for v in value]
        if isinstance(value,str) and value.startswith(str(OLD)+'/'):return str(OUT)+value[len(str(OLD)):]
        return value
    members=[relocate(m) for m in old['members']]
    for m in members:
        if m['harmonize_command'][2]!=str(P/'scripts/core_bounded_harmonizer_v3.py'):raise RuntimeError('EXACT_OLD_HELPER_REQUIRED')
        m['harmonize_command'][2]=str(P/'scripts/core_bounded_harmonizer_v4.py')
        if m['trait_id']=='sleep_apnea':
            at=m['harmonize_command'].index('--')
            m['harmonize_command'][at:at]=['--prepared-columns',str(manifest),'--expected-columns-sha256',manifest_sha]
            m['prepared_column_checkpoint']=checkpoint
    for name in ['core_apnea_empty_frame_diagnosis_controls_v1.py','core_apnea_empty_frame_diagnosis_controls_receipt_v1.json','core_apnea_empty_frame_scientific_diagnosis_v1.md']:
        bind(P/'reviews'/name)
    dependencies=dict(old['dependencies_sha256'])
    for name in ['core_bounded_harmonizer_v4.py','core_large35_checkpoint_adoption_v2.py','core_large35_checkpoint_continuation_v2.py','cleanup_core_bounded_spool_checkpoint_v1.py','prepare_core_large35_checkpoint_continuation_v2.py']:
        file=P/'scripts'/name;dependencies[str(file)]=sha(regular(file))
    dependencies.update(metadata)
    plan=dict(old)
    plan.update(schema='actual_core_large35_checkpoint_continuation_plan_v2',prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        private_namespace=str(OUT),members=members,dependencies_sha256=dependencies,
        checkpoint_metadata_sha256={**old['checkpoint_metadata_sha256'],**metadata},apnea_column_checkpoint=checkpoint,
        preserved_failed_checkpoint_plan=str(old_path),preserved_failed_checkpoint_plan_sha256=OLD_SHA,
        preserved_failed_checkpoint_execution_receipt=str(failed_path),preserved_failed_checkpoint_execution_sha256=metadata[str(failed_path)],
        zero_row_QC_correction='ONLY_EMPTY_FRAME_DROP_COUNT_ZERO;ALL_NONEMPTY_ORIGINAL_PREDICATES_UNCHANGED',
        disposal_helper=str(P/'scripts/cleanup_core_bounded_spool_checkpoint_v1.py'),
        new_source_worker_count=27,expected_new_worker_count=135,source_checkpoint_reused_as_new_worker=False,
        source_revalidation_qualification='Unchanged85/86 require current-plan source receipts; one guarded source check repeated, no source reparse or typed-column duplication.',
        capacity_last_actual_observation=record['resource_final'],
        capacity_qualification='Last genuine failed-worker snapshot, not a fresh meter or full-source RSS certification; unchanged live16GiB/all-core300GiB guards mandatory.',
        scientific_adapter_qualification='Existing component4 and v3 evidence retained; additional zero-row-only correction and read-only apnea checkpoint/disposal delta require independent review/root admission.',
        source_body_reads_in_preparation=0,output_body_reads_in_preparation=0,typed_piece_full_census_in_preparation=False,
        actual_small_frame_indices_inspected=[388,389],workers_launched=0,execution_admitted=False)
    adoption_gate(plan,current_content_hashes=False)
    if OUT.exists() or OUT.is_symlink():raise RuntimeError('NEW_SCIENTIFIC_CONTINUATION_NAMESPACE_MUST_BE_UNUSED')
    for m in members:
        for folder in ['prefilter','harmonized','munged','receipts']:(OUT/m['trait_id']/folder).mkdir(parents=True,exist_ok=False)
    for folder in ['receipts_v4','logs_v4','tmp','cache']:(OUT/folder).mkdir(parents=True,exist_ok=False)
    path=OUT/'core_large35_checkpoint_continuation_plan_v2.json';write_new(path,plan)
    print(json.dumps({'plan':str(path),'plan_sha256':sha(path),'remaining':27,'adopted':8,'new_workers':135,'checkpoint_manifest_sha256':manifest_sha,'workers_launched':0,'body_reads':0},indent=2))

if __name__=='__main__':main()

"""Metadata-only actual13-source plan; no bodies, runtime census or workers."""
import copy
from datetime import datetime,timezone
import json
from pathlib import Path
import importlib.util
from extension_replay_common_v4 import sha,write_new
from terminal_commit_common_v2 import require_committed
from validation_acquisition_producer_profile_v1 import regular,completed_acquisition

P=Path(__file__).resolve().parents[1]
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OUT=SSD/'validation_pipeline_replay_v4'
BASE=SSD/'sensitivities/sensitivity_operational_plan_v4_3_1.json'
BASE_SHA='05c0781aea2c7fe8f897d891e67bb741122e778382f35c8913946cf7469676b4'
ACQUISITION_SHA='3be0806d08c650bb733fcfb6d96f3f695aad2667f03d10d89de6a80ee8b715ee'
ADMISSION_SHA='52e8e9e91036318cd52c4b7c2b45f82cec3351e1c895fc39320ca935c053d8b4'


def prepare():
    if OUT.exists() or OUT.is_symlink():raise RuntimeError('NEW_VALIDATION_NAMESPACE_MUST_BE_UNUSED')
    metadata={}
    def bind(path,expected=None):
        path=regular(path);actual=sha(path)
        if expected is not None and actual!=expected:raise RuntimeError('ACTUAL_METADATA_OR_CODE_CHANGED: '+str(path))
        metadata[str(path)]=actual
        return actual
    bind(BASE,BASE_SHA);base=json.loads(BASE.read_text())
    design_path=P/'statistical_validation/validation_raw_replay_design_v1.json'
    bind(design_path,'f5291abd8989e6444f9aabe9ce2b0dd15258412acaa850dac1c59e72d86c7403')
    design=json.loads(design_path.read_text())
    acquisition_path=P/'manifests/validation_raw_acquisition_plan_v4_4.json';bind(acquisition_path,ACQUISITION_SHA)
    acquisition=json.loads(acquisition_path.read_text())
    admission_path=P/'manifests/validation_raw_acquisition_admission_v4_4.json';bind(admission_path,ADMISSION_SHA)
    executor=Path(acquisition['executor_path']);bind(executor,'6fec56a75ad71d841379436a461ac63d0d8cbcbe0ac62cf4aa960ead6a783d88')
    primary_paths=[P/'logs/validation_raw_acquisition_family_receipt_v4_4.json',Path(acquisition['private_namespace'])/'validation_raw_acquisition_family_receipt_v4_4.json']
    primary={str(path):bind(path) for path in primary_paths}
    if len(set(primary.values()))!=1:raise RuntimeError('TWO_GENUINE_FULL13_PRODUCER_COPIES_REQUIRED')
    binding=dict(plan_sha256=ACQUISITION_SHA,admission_sha256=ADMISSION_SHA,executor_sha256=metadata[str(executor)])
    bind(acquisition['terminal_seal_path']);require_committed(acquisition['pending_path'],acquisition['terminal_seal_path'],binding,primary)
    acquired=json.loads(primary_paths[0].read_text())
    if acquired['completed_source_count']!=13 or len(acquired['source_receipts'])!=13:raise RuntimeError('ACTUAL_COMPLETED_FULL13_PRODUCER_REQUIRED')
    individual=[];pairs=[]
    spec=importlib.util.spec_from_file_location('_exact108_validation_preparer',executor)
    producer=importlib.util.module_from_spec(spec);spec.loader.exec_module(producer)
    for member,item in zip(acquisition['members'],acquired['source_receipts']):
        origin=acquisition['receipt_origins'][member['source_id']]
        producer.source_receipt_gate(acquisition,member,item,ACQUISITION_SHA,verify_body=False)
        for path in [origin['primary_path'],origin['mirror_path']]:bind(path,item['sha256'])
        bind(origin['header_path'])
        individual.extend([origin['primary_path'],origin['mirror_path']])
        pairs.append(dict(member=member,primary_path=origin['primary_path'],mirror_path=origin['mirror_path'],sha256=item['sha256']))
    contract=dict(pending_path=acquisition['pending_path'],seal_path=acquisition['terminal_seal_path'],binding=binding,
                  primary_receipt_sha256=primary,metadata_sha256=dict(metadata),individual_receipt_paths=individual,
                  source_receipt_pairs=pairs,receipt_origins=acquisition['receipt_origins'],
                  producer_plan_path=str(acquisition_path),producer_executor_path=str(executor))
    dependencies=dict(base['dependencies_sha256']);dependencies.update(design['metadata_dependencies_sha256']);dependencies.update(acquisition['bound_sources'])
    for seal_name,expected in [('validation_collector_adapter_review_seal_v2.json','8938220336f9cb35a3e3396bff459c05093cf72add56dc3c8a8d757eb1a1753f'),('independent_original13_pipeline_prelaunch_seal_v3.json','82304e78550a6e277161a06ee470498556063a9e141baff992606c2039060be1')]:
        seal=P/'reviews'/seal_name;bind(seal,expected);r=json.loads(seal.read_text())
        for value in r.values():
            if isinstance(value,dict) and value and all(isinstance(x,str) and len(x)==64 for x in value.values()):
                for path,digest in value.items():
                    if path in dependencies and dependencies[path]!=digest:raise RuntimeError('CONTRADICTORY_INHERITED_REVIEW_IDENTITY')
                    dependencies[path]=digest
    # References, original processed streams and runtime bytes are inherited;
    # the unchanged adapter/validator/controller check them at actual execution.
    reference=design['reference_proposal']['reference_path'];alleles=design['reference_proposal']['allele_list_path']
    dependencies[reference]=design['reference_proposal']['reference_sha256_from_all13_original_receipts']
    dependencies[alleles]=design['reference_proposal']['allele_sha256_from_all13_original_receipts']
    collector=P.parent/'discovery_extension/scripts/47_stream_replication_sources.py'
    stream=P.parent/'discovery_extension/scripts/streaming_io.py';queue=P.parent/'discovery_extension/results/replication/replication_source_queue.tsv'
    for path in [collector,stream,queue,P/'scripts/96_replay_original_validation_collector_v2.py',P/'scripts/97_validate_original_validation_replay_v2.py',P/'scripts/101_run_original_validation_pipeline_v3.py',P/'scripts/102_prepare_original_validation_pipeline_v3.py']:
        bind(path,dependencies.get(str(path)))
    for name in ['validation_acquisition_producer_profile_v1.py','validation_collector_producer_profile_v1.py','run_original_validation_pipeline_v4.py','prepare_original_validation_pipeline_v4.py']:
        bind(P/'scripts'/name)
    python=design['collector_runtime']['executable'];logical=Path(python);runtime=logical.resolve(strict=True)
    runtime_sha=dependencies.pop(python)
    dependencies[str(runtime)]=runtime_sha
    profile=dict(logical_executable=python,resolved_binary=str(runtime),resolved_binary_sha256=runtime_sha,
                 declared_symlinks=[dict(path=str(v),target=str(v.readlink())) for v in [logical,*logical.parents] if v.is_symlink()],
                 historical_baseline_logical_binary_sha256=runtime_sha,historical_per_trait_binary_identity_claimed=False,
                 environment_read_only=True,current_runtime_bytes_freshly_hashed_in_preparation=False)
    members=[];workspace=OUT/'workspace'
    for original,raw in zip(design['members'],acquisition['members']):
        if original['source_id']!=raw['source_id'] or original['historical_source_body']['observed_sha256']!=raw['expected_sha256']:
            raise RuntimeError('EXACT_ORIGINAL13_SOURCE_IDENTITY_AND_ORDER_REQUIRED')
        sid=original['source_id'];row=original['frozen_source_metadata'];identity=dict(bytes=raw['expected_bytes'],md5=raw['expected_md5'],sha256=raw['expected_sha256'])
        for key in ['replication_munged_path','replication_receipt_path']:
            rel=Path(row[key])
            if rel.is_absolute() or '..' in rel.parts:raise RuntimeError('EXACT_ORIGINAL_RELATIVE_OUTPUT_ROUTES_REQUIRED')
        members.append(dict(source_id=sid,index=original['index'],body_path=raw['body_path'],original_source_identity=identity,
            frozen_source_metadata=row,new_munged=str(workspace/row['replication_munged_path']),new_receipt=str(workspace/row['replication_receipt_path']),
            new_qc=str(workspace/('discovery_extension/results/replication/qc/'+sid+'.tsv')),
            adapter_receipt=str(OUT/'receipts_v4'/(sid+'.adapter.json')),source_gate_receipt=str(OUT/'receipts_v4'/(sid+'.source.json')),
            comparison_receipt=str(OUT/'receipts_v4'/(sid+'.comparison.json')),original_munged=original['archived_processed_evidence']['path'],
            original_munged_sha256=original['historical_output_sha256'],original_qc=original['original_qc_path'],original_qc_sha256=original['original_qc_sha256'],
            original_receipt=original['original_receipt_path'],original_receipt_sha256=original['original_receipt_sha256'],
            expected_output_rows=int(original['original_qc']['output_rows']),effective_N_serialized_12g=original['effective_N_serialized_12g'],
            collector_argv=[str(collector),'--source-id',sid,'--queue',str(queue),'--reference',reference,'--hm3-alleles',alleles,'--acknowledge-network-gib',format(identity['bytes']/(1<<30),'.17g')]))
    dependencies.update(metadata)
    plan=copy.deepcopy(base)
    plan.update(schema='frozen_original13_validation_pipeline_plan_v4',prepared_utc=datetime.now(timezone.utc).isoformat(),members=members,member_count=13,jobs=[],
        new_estimator_fits=0,new_network_transfers=0,command_count=39,scope='ORIGINAL13_RAW_FILTER_SERIALIZATION_QC_REPLAY_ONLY',
        workspace=str(workspace),output_namespace=str(OUT),python=python,collector_runtime_resolved=str(runtime),collector_runtime_profile=profile,
        historical_per_trait_runtime_binary_identity_claimed=False,acquisition_plan_sha256=ACQUISITION_SHA,acquisition_terminal=contract,
        original_collector=str(collector),original_streaming_io=str(stream),queue=str(queue),reference=reference,alleles=alleles,
        adapter_path=str(P/'scripts/validation_collector_producer_profile_v1.py'),dependencies_sha256=dependencies,
        executor_path=str(P/'scripts/run_original_validation_pipeline_v4.py'),pending_path=str(OUT/'validation_pipeline_pending_v4.json'),
        terminal_seal_path=str(OUT/'validation_pipeline_terminal_v4.json'),primary_receipt_path=str(OUT/'validation_pipeline_execution_receipt_v4.json'),
        raw_acquisition_byte_count=10058648185,pipeline_output_ceiling_bytes=2<<30,execution_admitted=False,
        original_217_classification_and_41_native_estimates_unchanged=True,independent_replication_established=False,
        source_body_reads_in_preparation=0,reference_body_reads_in_preparation=0,runtime_censuses_in_preparation=0,workers_launched=0,
        preparation_qualification='CURRENT_GENUINE_PRODUCER_METADATA_AND_CODE_ONLY; source/current reference/archived output/runtime hashes inherited explicitly, not freshly reverified. Mandatory unchanged actual-execution gates remain.',
        required_review_filenames=['independent_original13_pipeline_prelaunch_v4.md','independent_original13_pipeline_prelaunch_v4.json','independent_original13_pipeline_prelaunch_seal_v4.json'])
    plan['guard'].update(worker_count=1,BLAS_threads=1,internal_floor_bytes=3<<30,SSD_floor_bytes=5<<30,observed_aggregate_worker_RSS_limit_bytes=2<<30,
        new_output_limit_bytes=2<<30,global_reservation_bytes=300<<30,shared_heavy_worker_lock=str(SSD/'native_heavy_worker.lock'),
        deadline_seconds=96*3600,per_worker_deadline_seconds=7200,poll_seconds=2,checkpoint_poll_seconds=2)
    completed_acquisition(plan)
    for path,digest in metadata.items():
        if sha(regular(path))!=digest:raise RuntimeError('CURRENT_PREPARATION_METADATA_CHANGED')
    for name in ['workspace','logs_v4','receipts_v4','tmp','cache']:(OUT/name).mkdir(parents=True,exist_ok=False)
    target=OUT/'validation_pipeline_plan_v4.json';write_new(target,plan)
    print(json.dumps(dict(plan=str(target),sha256=sha(target),sources=13,commands=39,workers_launched=0,body_reads=0),indent=2))


if __name__=='__main__':prepare()

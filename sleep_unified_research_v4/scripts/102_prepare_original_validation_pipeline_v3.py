#!/usr/bin/env python3
"""Freeze the original thirteen-source replay after actual acquisition commits."""
import copy
import hashlib
import json
from pathlib import Path
import subprocess
from terminal_commit_common_v2 import require_committed
from canonical_calibration_common_v4_5 import sha,utc,write_new

P=Path(__file__).resolve().parents[1]
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OUT=SSD/'validation_pipeline_replay_v3'
BASE=SSD/'sensitivities/sensitivity_operational_plan_v4_3_1.json'
BASE_SHA='05c0781aea2c7fe8f897d891e67bb741122e778382f35c8913946cf7469676b4'

def regular(path):
    p=Path(path)
    if p.is_symlink() or not p.is_file() or any(v.is_symlink() for v in p.parents):raise RuntimeError('REGULAR_PREPARATION_EVIDENCE_REQUIRED')
    return p

def prepare():
    if sha(BASE)!=BASE_SHA:raise RuntimeError('EXACT_NATIVE190_BASELINE_PLAN_REQUIRED')
    base=json.loads(BASE.read_text())
    design_path=P/'statistical_validation/validation_raw_replay_design_v1.json'
    if sha(regular(design_path))!='f5291abd8989e6444f9aabe9ce2b0dd15258412acaa850dac1c59e72d86c7403':raise RuntimeError('EXACT_ORIGINAL_VALIDATION_DESIGN_REQUIRED')
    design=json.loads(design_path.read_text())
    acquisition_path=P/'manifests/validation_raw_acquisition_plan_v4_3.json'
    if sha(regular(acquisition_path))!='8e93d548163bda284b6292b9780962060bbc3ff4afde44e5a916c80f3b69b708':raise RuntimeError('EXACT_ADMITTED_VALIDATION_ACQUISITION_PLAN_REQUIRED')
    acquisition=json.loads(acquisition_path.read_text())
    admission_path=P/'manifests/validation_raw_acquisition_admission_v4_3.json'
    if sha(regular(admission_path))!='812775f55aa2a2ee6f067c0ed5474f5592113732b99f375596312126724c067d':raise RuntimeError('EXACT_ACQUISITION_ADMISSION_REQUIRED')
    primary_paths=[P/'logs/validation_raw_acquisition_family_receipt_v4_3.json',Path(acquisition['private_namespace'])/'validation_raw_acquisition_family_receipt_v4_3.json']
    primary={str(path):sha(regular(path)) for path in primary_paths}
    if len(set(primary.values()))!=1:raise RuntimeError('TWO_EXACT_ACQUISITION_PRIMARY_COPIES_REQUIRED')
    binding=dict(plan_sha256=sha(acquisition_path),admission_sha256=sha(admission_path),executor_sha256=sha(acquisition['executor_path']))
    require_committed(acquisition['pending_path'],acquisition['terminal_seal_path'],binding,primary)
    acquired=json.loads(primary_paths[0].read_text())
    if acquired['status']!='ALL13_EXACT_ORIGINAL_SOURCE_BODIES_ACQUIRED' or acquired['completed_source_count']!=13 or len(acquired['source_receipts'])!=13:
        raise RuntimeError('ACTUAL_FULL13_ACQUISITION_PRODUCER_REQUIRED')
    dependencies=dict(base['dependencies_sha256']);dependencies.update(design['metadata_dependencies_sha256']);dependencies.update(acquisition['bound_sources'])
    metadata=[design_path,BASE,acquisition_path,admission_path,Path(acquisition['executor_path']),Path(acquisition['terminal_seal_path']),*primary_paths]
    individual=[];pairs=[]
    for member,item in zip(acquisition['members'],acquired['source_receipts']):
        path=regular(item['path']);mirror=regular(P/'source_provenance/validation_raw_acquisition_v4_3'/path.name)
        if sha(path)!=item['sha256'] or sha(mirror)!=item['sha256']:raise RuntimeError('EXACT_DOUBLE_SOURCE_RECEIPT_REQUIRED')
        r=json.loads(path.read_text())
        if r['member']!=member or r['status']!='EXACT_ORIGINAL_SOURCE_BODY_ACQUIRED' or r['plan_sha256']!=sha(acquisition_path) or not r['teardown']['teardown_verified'] or r['teardown']['remaining_group_members']:
            raise RuntimeError('EXACT_COMPLETED_SOURCE_PRODUCER_REQUIRED')
        individual.extend([str(path),str(mirror)]);metadata.extend([path,mirror])
        pairs.append(dict(member=member,primary_path=str(path),mirror_path=str(mirror),sha256=item['sha256']))
    for path in metadata:dependencies[str(path)]=sha(regular(path))
    contract=dict(pending_path=acquisition['pending_path'],seal_path=acquisition['terminal_seal_path'],binding=binding,
        primary_receipt_sha256=primary,metadata_sha256={str(path):sha(regular(path)) for path in metadata},individual_receipt_paths=individual,source_receipt_pairs=pairs)
    reference=design['reference_proposal']['reference_path'];alleles=design['reference_proposal']['allele_list_path']
    for path,digest in [(reference,design['reference_proposal']['reference_sha256_from_all13_original_receipts']),(alleles,design['reference_proposal']['allele_sha256_from_all13_original_receipts'])]:
        if sha(regular(path))!=digest:raise RuntimeError('FULL_ORIGINAL_COMMON_REFERENCE_IDENTITY_CHANGED')
        dependencies[path]=digest
    collector=P.parent/'discovery_extension/scripts/47_stream_replication_sources.py'
    stream=P.parent/'discovery_extension/scripts/streaming_io.py'
    queue=P.parent/'discovery_extension/results/replication/replication_source_queue.tsv'
    for path in [Path(__file__),P/'scripts/96_replay_original_validation_collector_v2.py',P/'scripts/101_run_original_validation_pipeline_v3.py',P/'scripts/97_validate_original_validation_replay_v2.py',P/'scripts/sensitivity_executor_v4_4.py',P/'scripts/terminal_commit_common_v2.py',P/'scripts/30_prepare_and_run_ssd_native_campaign.py',P/'scripts/native_stage_completion_v4_3.py',P/'scripts/canonical_calibration_common_v4_5.py',P/'scripts/canonical_calibration_common_v4_3.py',P/'scripts/extension_replay_common_v4.py',collector,stream,queue,P/'manifests/global_SSD_resource_reservation_v4_4.json',P/'logs/native190_root_independent_completion_addendum_v4.json']:
        dependencies[str(path)]=sha(regular(path))
    for name in ['validation13_original_collector_route_review_v1.md','validation13_original_collector_route_review_v1.json','validation13_original_collector_route_review_seal_v1.json']:
        path=P/'reviews'/name;dependencies[str(path)]=sha(regular(path))
    for name in ['92_replay_original_validation_collector_v1.py','93_prepare_original_validation_pipeline_v1.py','94_run_original_validation_pipeline_v1.py','95_validate_original_validation_replay_v1.py']:
        path=P/'scripts'/name;dependencies[str(path)]=sha(regular(path))
    review=P/'reviews/validation_collector_adapter_review_seal_v1.json'
    if sha(regular(review))!='729ed50869dc4e3f94b1032031d712974516e5d7b908054d47d2566d0395bef4':raise RuntimeError('PRESERVED_REJECTED_ADAPTER_REVIEW_REQUIRED')
    rejected=json.loads(review.read_text())
    maps=[v for k,v in rejected.items() if isinstance(v,dict) and v and all(isinstance(x,str) and len(x)==64 for x in v.values())]
    if len(maps)!=1:raise RuntimeError('EXACT_REJECTED_REVIEW_ARTIFACT_MAP_REQUIRED')
    for path,digest in maps[0].items():
        if sha(regular(path))!=digest:raise RuntimeError('REJECTED_REVIEW_ARTIFACT_CHANGED')
        dependencies[path]=digest
    dependencies[str(review)]=sha(review)
    # Preserve the rejected whole-pipeline oracle and independently accepted science component.
    for oldname,oldsha in [('98_prepare_original_validation_pipeline_v2.py','eb80c1e2185490c921bb1e87d0e5ca45fe4cc788664e9042e90caec845b87a7c'),('99_run_original_validation_pipeline_v2.py','c47dccbaf3e353c28a5cec4b51e9a191b4d828c93d29c3c4837cb0f327eee416')]:
        oldpath=P/'scripts'/oldname
        if sha(regular(oldpath))!=oldsha:raise RuntimeError('PRESERVED_REJECTED_PIPELINE_CODE_CHANGED')
        dependencies[str(oldpath)]=oldsha
    for sealname,sealdigest in [('independent_original13_pipeline_prelaunch_seal_v2.json','dc3c560ada5086c778d51af4063091447846ad7269963ac6efbe2996a5cc4730'),('validation_collector_adapter_review_seal_v2.json','8938220336f9cb35a3e3396bff459c05093cf72add56dc3c8a8d757eb1a1753f')]:
        sealpath=P/'reviews'/sealname
        if sha(regular(sealpath))!=sealdigest:raise RuntimeError('EXACT_PRESERVED_REVIEW_SEAL_REQUIRED')
        sealed=json.loads(sealpath.read_text())
        for mapping in sealed.values():
            if isinstance(mapping,dict) and mapping and all(isinstance(value,str) and len(value)==64 for value in mapping.values()):
                for path,digest in mapping.items():
                    if sha(regular(path))!=digest:raise RuntimeError('PRESERVED_REVIEW_ARTIFACT_CHANGED')
                    if path in dependencies and dependencies[path]!=digest:raise RuntimeError('CONTRADICTORY_REVIEW_DEPENDENCY')
                    dependencies[path]=digest
        if sha(regular(sealpath))!=sealdigest:raise RuntimeError('REVIEW_SEAL_CHANGED_DURING_CONSUMPTION')
        dependencies[str(sealpath)]=sealdigest
    python=design['collector_runtime']['executable']
    if subprocess.run([python,'-B','--version'],capture_output=True,text=True,check=True).stdout.strip()!='Python 3.9.23':raise RuntimeError('DECLARED_ORIGINAL_COLLECTOR_RUNTIME_CHANGED')
    runtime=Path(python).resolve(strict=True);runtime_sha=sha(regular(runtime))
    if dependencies.get(python)!=runtime_sha:raise RuntimeError('BASELINE_DECLARED_RUNTIME_BINARY_IDENTITY_REQUIRED')
    baseline_runtime_sha=dependencies.pop(python)
    dependencies[str(runtime)]=runtime_sha
    logical=Path(python)
    runtime_profile=dict(logical_executable=python,resolved_binary=str(runtime),resolved_binary_sha256=runtime_sha,
        declared_symlinks=[{'path':str(v),'target':str(v.readlink())} for v in [logical,*logical.parents] if v.is_symlink()],
        historical_baseline_logical_binary_sha256=baseline_runtime_sha,
        historical_per_trait_binary_identity_claimed=False,environment_read_only=True)
    # The known interpreter link is explicit; all scientific input dependencies must be regular.
    for path,digest in dependencies.items():
        if sha(regular(path))!=digest:raise RuntimeError('REGULAR_FULL_PREPARATION_DEPENDENCY_REQUIRED')
    members=[];workspace=OUT/'workspace'
    for original,raw in zip(design['members'],acquisition['members']):
        if original['source_id']!=raw['source_id'] or original['historical_source_body']['observed_sha256']!=raw['expected_sha256']:raise RuntimeError('UNCHANGED_SOURCE_ORDER_AND_IDENTITY_REQUIRED')
        sid=original['source_id'];row=original['frozen_source_metadata'];identity=dict(bytes=raw['expected_bytes'],md5=raw['expected_md5'],sha256=raw['expected_sha256'])
        body=regular(raw['body_path'])
        md5=hashlib.md5();digest=hashlib.sha256()
        with body.open('rb') as f:
            for b in iter(lambda:f.read(65536),b''):md5.update(b);digest.update(b)
        if body.stat().st_size!=identity['bytes'] or (md5.hexdigest(),digest.hexdigest())!=(identity['md5'],identity['sha256']):raise RuntimeError('ACTUAL_CURRENT_FULL13_BODIES_REQUIRED')
        archive=original['archived_processed_evidence']['path']
        if sha(regular(archive))!=original['historical_output_sha256']:raise RuntimeError('CURRENT_ORIGINAL_PROCESSED_IDENTITY_REQUIRED')
        for key in ['replication_munged_path','replication_receipt_path']:
            rel=Path(row[key])
            if rel.is_absolute() or '..' in rel.parts:raise RuntimeError('ORIGINAL_RELATIVE_PRIVATE_ROUTE_REQUIRED')
        member=dict(source_id=sid,index=original['index'],body_path=raw['body_path'],original_source_identity=identity,
            frozen_source_metadata=row,new_munged=str(workspace/row['replication_munged_path']),
            new_receipt=str(workspace/row['replication_receipt_path']),new_qc=str(workspace/('discovery_extension/results/replication/qc/'+sid+'.tsv')),
            adapter_receipt=str(OUT/'receipts_v4'/(sid+'.adapter.json')),source_gate_receipt=str(OUT/'receipts_v4'/(sid+'.source.json')),
            comparison_receipt=str(OUT/'receipts_v4'/(sid+'.comparison.json')),original_munged=archive,
            original_munged_sha256=original['historical_output_sha256'],original_qc=original['original_qc_path'],original_qc_sha256=original['original_qc_sha256'],
            original_receipt=original['original_receipt_path'],original_receipt_sha256=original['original_receipt_sha256'],
            expected_output_rows=int(original['original_qc']['output_rows']),effective_N_serialized_12g=original['effective_N_serialized_12g'],
            collector_argv=[str(collector),'--source-id',sid,'--queue',str(queue),'--reference',reference,'--hm3-alleles',alleles,'--acknowledge-network-gib',format(identity['bytes']/(1<<30),'.17g')])
        members.append(member)
    plan=copy.deepcopy(base)
    plan.update(schema='frozen_original13_validation_pipeline_plan_v3',prepared_utc=utc(),members=members,member_count=13,jobs=[],
        new_estimator_fits=0,new_network_transfers=0,command_count=39,scope='ORIGINAL13_RAW_FILTER_SERIALIZATION_QC_REPLAY_ONLY',
        workspace=str(workspace),output_namespace=str(OUT),python=python,collector_runtime_resolved=str(runtime),collector_runtime_profile=runtime_profile,
        historical_per_trait_runtime_binary_identity_claimed=False,acquisition_plan_sha256=sha(acquisition_path),acquisition_terminal=contract,
        original_collector=str(collector),original_streaming_io=str(stream),queue=str(queue),reference=reference,alleles=alleles,
        adapter_path=str(P/'scripts/96_replay_original_validation_collector_v2.py'),dependencies_sha256=dependencies,
        executor_path=str(P/'scripts/101_run_original_validation_pipeline_v3.py'),
        pending_path=str(OUT/'validation_pipeline_pending_v3.json'),terminal_seal_path=str(OUT/'validation_pipeline_terminal_v3.json'),
        primary_receipt_path=str(OUT/'validation_pipeline_execution_receipt_v3.json'),
        raw_acquisition_byte_count=10058648185,pipeline_output_ceiling_bytes=2<<30,
        original_217_classification_and_41_native_estimates_unchanged=True,independent_replication_established=False,
        required_review_filenames=['independent_original13_pipeline_prelaunch_v3.md','independent_original13_pipeline_prelaunch_v3.json','independent_original13_pipeline_prelaunch_seal_v3.json'])
    plan['guard'].update(worker_count=1,BLAS_threads=1,internal_floor_bytes=3<<30,SSD_floor_bytes=5<<30,
        observed_aggregate_worker_RSS_limit_bytes=2<<30,new_output_limit_bytes=2<<30,global_reservation_bytes=300<<30,
        shared_heavy_worker_lock=str(SSD/'native_heavy_worker.lock'),deadline_seconds=96*3600,per_worker_deadline_seconds=7200,poll_seconds=2,checkpoint_poll_seconds=2)
    for path,digest in dependencies.items():
        if sha(regular(path))!=digest:raise RuntimeError('PREPARATION_DEPENDENCY_CHANGED')
    import importlib.util
    spec=importlib.util.spec_from_file_location('_validation_preparation_consumer',plan['adapter_path'])
    adapter=importlib.util.module_from_spec(spec);spec.loader.exec_module(adapter)
    adapter.runtime_gate(plan);adapter.completed_acquisition(plan)
    require_committed(contract['pending_path'],contract['seal_path'],contract['binding'],primary)
    for directory in ['workspace','logs_v4','receipts_v4','tmp','cache']:(OUT/directory).mkdir(parents=True,exist_ok=True)
    target=OUT/'validation_pipeline_plan_v3.json';write_new(target,plan)
    print(json.dumps(dict(plan=str(target),sha256=sha(target),sources=13,commands=39,workers_launched=0)))

if __name__=='__main__':prepare()

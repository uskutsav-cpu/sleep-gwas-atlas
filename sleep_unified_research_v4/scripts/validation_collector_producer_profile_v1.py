#!/usr/bin/env python3
"""Run unchanged original validation filtering on one exact local source body."""
import argparse
from contextlib import contextmanager
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from terminal_commit_common_v2 import require_committed

def hashes(path):
    md5=hashlib.md5();sha=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(65536),b''):md5.update(block);sha.update(block)
    return md5.hexdigest(),sha.hexdigest()

def regular(path):
    p=Path(path)
    if p.is_symlink() or not p.is_file() or any(parent.is_symlink() for parent in p.parents):
        raise RuntimeError('REGULAR_FROZEN_VALIDATION_INPUT_REQUIRED: '+str(p))
    return p

def fixed_hashes(mapping):
    if not mapping:raise RuntimeError('NONEMPTY_FROZEN_VALIDATION_METADATA_REQUIRED')
    for path,digest in mapping.items():
        if hashes(regular(path))[1]!=digest:raise RuntimeError('FROZEN_VALIDATION_INPUT_CHANGED: '+path)

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def runtime_gate(plan):
    profile=plan['collector_runtime_profile']
    logical=Path(plan['python']);resolved=regular(profile['resolved_binary'])
    if str(logical)!=profile['logical_executable'] or str(logical.resolve(strict=True))!=str(resolved):
        raise RuntimeError('EXACT_DECLARED_COLLECTOR_RUNTIME_ROUTE_REQUIRED')
    actual=[{'path':str(v),'target':str(v.readlink())} for v in [logical,*logical.parents] if v.is_symlink()]
    if actual!=profile['declared_symlinks'] or hashes(resolved)[1]!=profile['resolved_binary_sha256'] or plan['dependencies_sha256'].get(str(resolved))!=profile['resolved_binary_sha256']:
        raise RuntimeError('DECLARED_READ_ONLY_RUNTIME_IDENTITY_CHANGED')

from validation_acquisition_producer_profile_v1 import completed_acquisition

def run(plan_path,expected_plan_sha,source_id):
    if hashes(regular(plan_path))[1]!=expected_plan_sha:raise RuntimeError('EXACT_FROZEN_VALIDATION_REPLAY_PLAN_REQUIRED')
    plan=json.loads(Path(plan_path).read_text())
    if plan['member_count']!=13 or len(plan['members'])!=13 or plan['new_estimator_fits']!=0:
        raise RuntimeError('UNCHANGED_THIRTEEN_SOURCE_PREPROCESSING_SCOPE_REQUIRED')
    selected=[m for m in plan['members'] if m['source_id']==source_id]
    if len(selected)!=1:raise RuntimeError('EXACTLY_ONE_FROZEN_VALIDATION_SOURCE_REQUIRED')
    member=selected[0];fixed_hashes(plan['dependencies_sha256']);runtime_gate(plan)
    completed_acquisition(plan)
    body=regular(member['body_path'])
    original=member['original_source_identity']
    def current_body():
        if body.stat().st_size!=original['bytes'] or hashes(regular(body))!=(original['md5'],original['sha256']):
            raise RuntimeError('FULL_ORIGINAL_VALIDATION_BODY_CHANGED')
    workspace=Path(plan['workspace'])
    if not workspace.is_absolute() or not workspace.is_dir() or workspace.is_symlink() or any(p.is_symlink() for p in workspace.parents):
        raise RuntimeError('REGULAR_PRIVATE_VALIDATION_WORKSPACE_REQUIRED')
    row=member['frozen_source_metadata']
    for key,actual_key in [('replication_munged_path','new_munged'),('replication_receipt_path','new_receipt')]:
        relative=Path(row[key])
        if relative.is_absolute() or '..' in relative.parts or str(workspace/relative)!=member[actual_key]:
            raise RuntimeError('EXACT_PRIVATE_ORIGINAL_OUTPUT_ROUTE_REQUIRED')
    expected_qc=workspace/('discovery_extension/results/replication/qc/'+source_id+'.tsv')
    if member['new_qc']!=str(expected_qc):raise RuntimeError('EXACT_PRIVATE_ORIGINAL_QC_ROUTE_REQUIRED')
    namespace=Path(plan['output_namespace'])
    if str(workspace)!=str(namespace/'workspace') or namespace.is_symlink() or any(v.is_symlink() for v in namespace.parents):
        raise RuntimeError('EXACT_PRIVATE_WORKSPACE_NAMESPACE_REQUIRED')
    if member['adapter_receipt']!=str(namespace/'receipts_v4'/(source_id+'.adapter.json')):
        raise RuntimeError('EXACT_PRIVATE_ADAPTER_RECEIPT_ROUTE_REQUIRED')
    for path in [member['new_munged'],member['new_receipt'],member['new_qc'],member['adapter_receipt']]:
        target=Path(path)
        if target.exists() or target.is_symlink() or any(parent.is_symlink() for parent in target.parents):
            raise RuntimeError('PRIOR_PRIVATE_VALIDATION_OUTPUT_PRESERVED_NO_RETRY')
    current_body()
    collector=regular(plan['original_collector']);stream_code=regular(plan['original_streaming_io'])
    stream=load('_original_validation_streaming',stream_code)
    # Original import resolves to this exact unmodified module, with no network call.
    prior_stream=sys.modules.get('streaming_io');sys.modules['streaming_io']=stream
    try:module=load('_original_validation_collector',collector)
    finally:
        if prior_stream is None:sys.modules.pop('streaming_io',None)
        else:sys.modules['streaming_io']=prior_stream
    reference=regular(plan['reference']);alleles=regular(plan['alleles']);queue=regular(plan['queue'])
    routes={'discovery_extension/scripts/streaming_io.py':stream_code,
            'discovery_extension/scripts/47_stream_replication_sources.py':collector,
            'discovery_extension/results/replication/replication_source_queue.tsv':queue,
            'discovery_extension/data/reference/panukbb_hm3_variant_reference.tsv.gz':reference,
            'ref/w_hm3.snplist':alleles}
    def scoped_Path(*args):
        path=Path(*args)
        if path.is_absolute():
            temporary=(path.parent==Path(member['new_munged']).parent and path.name.startswith(Path(member['new_munged']).name+'.') and path.name.endswith('.tmp'))
            if path not in [collector,stream_code,reference,alleles,queue] and not temporary:raise RuntimeError('UNREGISTERED_ABSOLUTE_COLLECTOR_PATH')
            if '..' in path.parts or any(parent.is_symlink() for parent in path.parents):raise RuntimeError('COLLECTOR_PATH_PARENT_ESCAPE')
            return path
        if str(path) not in routes:raise RuntimeError('UNREGISTERED_RELATIVE_COLLECTOR_PATH: '+str(path))
        return routes[str(path)]
    def frozen_head(url):
        if url!=row['replication_source_url']:raise RuntimeError('UNREGISTERED_COLLECTOR_SOURCE_URL')
        # This is frozen queue metadata, not an assertion about current HTTP headers.
        # Actual transport is separately bound by the completed acquisition producer.
        fields={'content-length':row['replication_content_length_bytes'],'etag':row['replication_etag']}
        if row['replication_storage_mode']=='VERSIONED_REMOTE_STREAMING':fields['x-goog-generation']=row['replication_source_generation']
        return fields
    @contextmanager
    def local_stream(*,url,expected_md5,expected_size_bytes,timeout_seconds):
        if url!=row['replication_source_url'] or expected_md5!=row['replication_checksum'] or expected_size_bytes!=original['bytes'] or timeout_seconds!=240:
            raise RuntimeError('UNREGISTERED_ORIGINAL_STREAM_ARGUMENTS')
        with stream.open_verified_gzip_text(local_path=body,expected_md5=expected_md5,
                 expected_size_bytes=expected_size_bytes,timeout_seconds=timeout_seconds) as (text,receipt):
            receipt.update(locked_original_remote_url=url,replay_route='FULL_SHA_VERIFIED_PRIVATE_LOCAL_BODY',
                local_transport_does_not_claim_original_HTTP_headers=True,
                acquisition_plan_sha256=plan['acquisition_plan_sha256'])
            yield text,receipt
        if receipt['observed_sha256']!=original['sha256'] or receipt['verification_status']!='PASS':
            raise RuntimeError('FULL_LOCAL_GZIP_EOF_IDENTITY_CHANGED')
    expected_argv=[str(collector),'--source-id',source_id,'--queue',str(queue),'--reference',str(reference),
        '--hm3-alleles',str(alleles),'--acknowledge-network-gib',format(original['bytes']/(1<<30),'.17g')]
    if member['collector_argv']!=expected_argv:raise RuntimeError('FROZEN_COLLECTOR_ARGUMENTS_CHANGED')
    module.ROOT=workspace;module.Path=scoped_Path;module.head=frozen_head;module.open_verified_gzip_text=local_stream
    old_argv=sys.argv;sys.argv=expected_argv
    try:module.main()
    finally:sys.argv=old_argv
    current_body();fixed_hashes(plan['dependencies_sha256'])
    completed_acquisition(plan)
    if hashes(regular(plan_path))[1]!=expected_plan_sha:raise RuntimeError('VALIDATION_REPLAY_PLAN_CHANGED_AFTER_COLLECTOR')
    consumed_outputs={path:hashes(regular(path))[1] for path in [member['new_munged'],member['new_receipt'],member['new_qc']]}
    fixed_hashes(consumed_outputs)
    result=json.loads(Path(member['new_receipt']).read_text())
    fixed_hashes(consumed_outputs)
    if result['pipeline_status']!='STREAM_HARMONIZE_DIRECT_MUNGE_PASS' or result['replication_source_id']!=source_id or result['source_url']!=row['replication_source_url'] or result['source_verification']['observed_sha256']!=original['sha256']:
        raise RuntimeError('EXACT_ORIGINAL_COLLECTOR_RESULT_REQUIRED')
    sidecar=dict(schema='original_validation_collector_local_replay_sidecar_v2',plan_sha256=expected_plan_sha,
        source_id=source_id,adapter_sha256=hashes(Path(__file__))[1],original_source_identity=original,
        full_local_source_retained=True,original_URL_and_queue_preserved=True,
        HTTP_HEAD_not_requeried_in_collector=True,offline_gate='Completed separately admitted original-body acquisition; original queue metadata supplies the legacy collector HEAD check.',
        legacy_collector_receipt_fields_preserved='Original full_resolution_local_retention:false and remote-identity text are historical collector literals; the acquired source is retained locally by this replay.',
        full_source_gzip_CRC_EOF_MD5_size_SHA256_verified=True,original_filters_and_serialization_unchanged=True,
        output_sha256=consumed_outputs,
        source_or_independent_replication_admitted=False,new_estimator_fits=0)
    payload=json.dumps(sidecar,indent=2,allow_nan=False)+'\n';expected=hashlib.sha256(payload.encode()).hexdigest()
    fixed_hashes(consumed_outputs);runtime_gate(plan)
    if hashes(regular(plan_path))[1]!=expected_plan_sha:raise RuntimeError('FIXED_PLAN_REQUIRED_BEFORE_SIDECAR')
    with Path(member['adapter_receipt']).open('x') as f:f.write(payload)
    if hashes(regular(member['adapter_receipt']))[1]!=expected:raise RuntimeError('INTENDED_ADAPTER_RECEIPT_CHANGED')
    current_body();completed_acquisition(plan);fixed_hashes(plan['dependencies_sha256']);runtime_gate(plan)
    fixed_hashes(consumed_outputs)
    if hashes(regular(plan_path))[1]!=expected_plan_sha:raise RuntimeError('FIXED_PLAN_REQUIRED_AFTER_SIDECAR')
    if hashes(regular(member['adapter_receipt']))[1]!=expected:raise RuntimeError('FIXED_INTENDED_SIDECAR_REQUIRED_AT_RETURN')

def main():
    p=argparse.ArgumentParser();p.add_argument('--plan',type=Path,required=True);p.add_argument('--plan-sha256',required=True);p.add_argument('--source-id',required=True)
    a=p.parse_args();run(a.plan,a.plan_sha256,a.source_id)

if __name__=='__main__':main()

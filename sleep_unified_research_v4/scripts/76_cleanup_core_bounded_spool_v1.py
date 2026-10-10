#!/usr/bin/env python3
"""Dispose only owned generated scratch after exact original-content/QC proof.

No source, original output, estimator or historical archive is altered. The
supervising controller must retain its shared lock and resource monitoring.
"""
import argparse
import datetime
import json
from pathlib import Path
import re

from extension_replay_common_v3 import sha, write_new, check_bindings

POLICY = 'DELETE_ONLY_OWN_GENERATED_SCRATCH_AFTER_EXACT_FULL_ORIGINAL_CONTENT_QC_PASS;FREEZE_COPIED_CANDIDATE_MANIFESTS_AND_ENTIRE_INVENTORY_BEFORE_DISPOSAL;PRESERVE_FAILED_SCRATCH_AND_ALL_ORIGINAL_INPUTS'


def regular(path):
    p = Path(path)
    if p.is_symlink() or not p.is_file() or any(x.is_symlink() for x in p.parents):
        raise RuntimeError('OWNED_SCRATCH_EVIDENCE_NOT_REGULAR: '+str(p))
    return p


def admissible(relative):
    # Include the SSD's transport sidecars only when their exact partner name
    # belongs to this generated spool. No unrelated hidden files are accepted.
    parts = list(relative.parts)
    if parts[-1].startswith('._'):
        parts[-1] = parts[-1][2:]
    name = '/'.join(parts)
    return name in {'global_harmonization_candidate_receipt.json',
                    'columns/column_preparation_manifest.json',
                    'required_map_keys.sqlite3', 'first_surviving_SNP.sqlite3'} or bool(
        re.fullmatch(r'columns/(?:typed|tokens)_[0-9]{2}_[0-9]{8}\.pkl', name)
        or re.fullmatch(r'columns/coordinate_tokens_(?:CHR|BP)_[0-9]{8}\.pkl', name)
        or re.fullmatch(r'(?:post_liftover|post_ordinary_QC)/[0-9]{8}\.pkl', name))


def run(plan_path, plan_hash, trait):
    if sha(plan_path) != plan_hash:
        raise RuntimeError('FROZEN_CORE_PLAN_CHANGED')
    plan = json.loads(plan_path.read_text())
    check_bindings(plan)
    if plan['ephemeral_spool_policy'] != POLICY:
        raise RuntimeError('EXACT_OWNED_SCRATCH_POLICY_REQUIRED')
    member = next(m for m in plan['members'] if m['trait_id'] == trait)
    spool = Path(member['ephemeral_spool'])
    expected = Path(plan['private_namespace'])/trait/'ephemeral_bounded_spool'
    if spool != expected or spool.is_symlink() or not spool.is_dir() or any(x.is_symlink() for x in spool.parents):
        raise RuntimeError('EXACT_PRIVATE_GENERATED_SPOOL_REQUIRED')
    candidate_path = regular(spool/'global_harmonization_candidate_receipt.json')
    column_path = regular(spool/'columns/column_preparation_manifest.json')
    candidate = json.loads(candidate_path.read_text())
    columns = json.loads(column_path.read_text())
    comparison_path = regular(member['comparison_receipt'])
    source_path = regular(member['source_gate_receipt'])
    comparison = json.loads(comparison_path.read_text())
    source = json.loads(source_path.read_text())
    raw = member['original_design']['raw']
    if (comparison['status'] != 'EXACT_CORE_HARMONIZED_MUNGED_CONTENT_AND_QC_REPLAY_PASS'
        or source['status'] != 'EXACT_CORE_RAW_SOURCE_GATE_PASS'
        or comparison['trait'] != trait or source['trait'] != trait
        or comparison['plan_sha256'] != plan_hash or source['plan_sha256'] != plan_hash
        or comparison['raw_sha256'] != raw['sealed_verified_sha256']
        or source['raw_sha256'] != raw['sealed_verified_sha256']):
        raise RuntimeError('EXACT_WHOLE_ORIGINAL_CONTENT_QC_AND_SOURCE_PROOF_REQUIRED_BEFORE_DISPOSAL')
    if candidate['status'] != 'CANDIDATE_GLOBAL_HARMONIZATION_PRODUCED_NOT_SCIENTIFICALLY_ADMITTED':
        raise RuntimeError('EXACT_GENERATED_CANDIDATE_RECEIPT_REQUIRED')
    if (candidate['source_sha256'] != raw['sealed_verified_sha256']
        or candidate['expected_rows'] != member['original_design']['harmonizer_input_rows']
        or candidate['rows_out'] != member['original_design']['harmonized_rows']
        or candidate['output_sha256'] != sha(regular(member['harmonized']))
        or candidate['column_manifest_sha256'] != sha(column_path)
        or columns['source_sha256_before'] != raw['sealed_verified_sha256']
        or columns['source_sha256_after'] != raw['sealed_verified_sha256']):
        raise RuntimeError('GENERATED_SPOOL_AND_EXACT_FINAL_OUTPUT_BINDINGS_DIFFER')
    preconditions = {str(p):sha(p) for p in [comparison_path, source_path,
        regular(member['harmonized']), regular(member['munged']),
        regular(member['harmonization_qc'])]}
    pieces = [item for mapping in [columns['typed_pieces'],
        columns['preinference_token_pieces'], columns.get('coordinate_token_pieces', {})]
        for items in mapping.values() for item in items]
    pieces += candidate['post_liftover_frames'] + candidate['post_ordinary_QC_frames']
    for item in pieces:
        path = regular(item['path'])
        if not path.is_relative_to(spool) or sha(path) != item['sha256']:
            raise RuntimeError('CONSUMED_GENERATED_PIECE_CHANGED_BEFORE_DISPOSAL')
    files = sorted(p for p in spool.rglob('*') if p.is_file() or p.is_symlink())
    if any(not admissible(p.relative_to(spool)) for p in files):
        raise RuntimeError('UNREGISTERED_FILE_PRESERVED_NO_SCRATCH_DISPOSAL')
    inventory = {str(regular(p)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in files}
    for path, value in preconditions.items():
        if sha(regular(path)) != value:
            raise RuntimeError('EXACT_CONTENT_QC_PROOF_CHANGED_BEFORE_DISPOSAL')
    if sha(plan_path) != plan_hash:
        raise RuntimeError('CORE_PLAN_CHANGED_BEFORE_DISPOSAL')
    check_bindings(plan)
    write_new(member['bounded_receipt'], candidate)
    write_new(member['column_manifest'], columns)
    frozen = {'schema':'owned_core_ephemeral_spool_inventory_v1',
        'plan_sha256':plan_hash,'trait':trait,'source_namespace':str(spool),
        'scientific_success_proof_sha256':preconditions,'file_inventory':inventory,
        'generated_bytes':sum(v['bytes'] for v in inventory.values()),
        'deletion_policy':POLICY,'copied_candidate_receipt_sha256':sha(member['bounded_receipt']),
        'copied_column_manifest_sha256':sha(member['column_manifest'])}
    write_new(member['spool_inventory'], frozen)
    # The complete immutable intent/inventory is durable before the first
    # unlink. Interruption leaves this evidence and remaining scratch intact.
    deleted = 0
    for path, item in inventory.items():
        file = regular(path)
        if file.stat().st_size != item['bytes'] or sha(file) != item['sha256']:
            raise RuntimeError('GENERATED_SCRATCH_CHANGED_DURING_DISPOSAL')
        file.unlink(); deleted += 1
    for folder in sorted((p for p in spool.rglob('*') if p.is_dir()), key=lambda p:len(p.parts), reverse=True):
        if folder.is_symlink():raise RuntimeError('SPOOL_DIRECTORY_SYMLINK_PRESERVED')
        folder.rmdir()
    spool.rmdir()
    for path, value in preconditions.items():
        if sha(regular(path)) != value:
            raise RuntimeError('EXACT_FINAL_OUTPUT_CHANGED_AFTER_DISPOSAL')
    if sha(plan_path) != plan_hash:raise RuntimeError('CORE_PLAN_CHANGED_AFTER_DISPOSAL')
    check_bindings(plan)
    record = {'schema':'owned_core_ephemeral_spool_cleanup_v1',
        'status':'OWN_GENERATED_SCRATCH_DISPOSED_AFTER_EXACT_FULL_CONTENT_QC_PASS',
        'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'plan_sha256':plan_hash,'trait':trait,'deleted_files':deleted,
        'deleted_bytes':frozen['generated_bytes'],'spool_absent':not spool.exists(),
        'inventory_sha256':sha(member['spool_inventory']),
        'scientific_success_proof_sha256':preconditions,
        'copied_candidate_receipt_sha256':sha(member['bounded_receipt']),
        'copied_column_manifest_sha256':sha(member['column_manifest']),
        'original_input_or_output_files_deleted':0}
    write_new(member['spool_cleanup_receipt'], record)
    print(json.dumps({'trait':trait,'status':record['status'],'generated_bytes_disposed':record['deleted_bytes']}), flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--plan-sha',required=True);p.add_argument('--trait',required=True)
    a=p.parse_args();run(a.plan,a.plan_sha,a.trait)


if __name__=='__main__':main()

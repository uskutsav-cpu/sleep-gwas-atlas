#!/usr/bin/env python3
"""Dispose the exact registered private spool after frozen whole-content proof."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import re

from extension_replay_common_v3 import sha, write_new, check_bindings
from terminal_commit_common_v2 import sync_directory

POLICY = 'DELETE_ONLY_OWN_GENERATED_SCRATCH_AFTER_EXACT_FULL_ORIGINAL_CONTENT_QC_PASS;FREEZE_COPIED_CANDIDATE_MANIFESTS_AND_ENTIRE_INVENTORY_BEFORE_DISPOSAL;PRESERVE_FAILED_SCRATCH_AND_ALL_ORIGINAL_INPUTS'
P = Path(__file__).resolve().parents[1]


def regular(path):
    p = Path(path)
    if p.is_symlink() or not p.is_file() or any(x.is_symlink() for x in p.parents):
        raise RuntimeError('EXACT_PRIVATE_EVIDENCE_NOT_REGULAR: ' + str(p))
    return p


def unchanged(mapping):
    for path, digest in mapping.items():
        if sha(regular(path)) != digest:
            raise RuntimeError('FROZEN_DISPOSAL_EVIDENCE_CHANGED: ' + path)


def copy_exact(source, destination, expected):
    source = regular(source)
    destination = Path(destination)
    if destination.exists() or destination.is_symlink() or any(x.is_symlink() for x in destination.parents):
        raise RuntimeError('PRIOR_OR_REDIRECTED_RETAINED_EVIDENCE_PRESERVED')
    if sha(source) != expected:
        raise RuntimeError('FROZEN_GENERATED_METADATA_CHANGED_BEFORE_COPY')
    with source.open('rb') as src, destination.open('xb') as dst:
        for block in iter(lambda: src.read(65536), b''):
            dst.write(block)
        dst.flush()
        os.fsync(dst.fileno())
    if sha(regular(destination)) != expected or sha(source) != expected:
        raise RuntimeError('EXACT_GENERATED_METADATA_COPY_DIFFERS')


def run(plan_path, plan_hash, trait, candidate_hash, columns_hash, comparison_worker_hash):
    plan_path = regular(plan_path)
    if sha(plan_path) != plan_hash:
        raise RuntimeError('FROZEN_CORE_PLAN_CHANGED')
    plan = json.loads(plan_path.read_text())
    check_bindings(plan)
    root = Path(plan['private_namespace'])
    if root != plan_path.parent or any(x.is_symlink() for x in [root, *root.parents]):
        raise RuntimeError('EXACT_PRIVATE_PLAN_NAMESPACE_REQUIRED')
    if not re.fullmatch(r'[A-Za-z0-9_]+', trait) or plan['ephemeral_spool_policy'] != POLICY:
        raise RuntimeError('EXACT_REGISTERED_TRAIT_AND_SCRATCH_POLICY_REQUIRED')
    matches = [m for m in plan['members'] if m['trait_id'] == trait]
    if len(matches) != 1:
        raise RuntimeError('EXACT_SINGLE_FROZEN_TRAIT_REQUIRED')
    member = matches[0]
    folder = root / trait
    spool = folder / 'ephemeral_bounded_spool'
    if Path(member['ephemeral_spool']) != spool or spool.is_symlink() or not spool.is_dir() or any(x.is_symlink() for x in spool.parents):
        raise RuntimeError('EXACT_PRIVATE_GENERATED_SPOOL_REQUIRED')
    retained_names = {'bounded_receipt': 'bounded_harmonization.json',
                      'column_manifest': 'native_column_preparation.json',
                      'spool_inventory': 'ephemeral_spool_inventory.json',
                      'spool_cleanup_receipt': 'ephemeral_spool_cleanup.json',
                      'source_gate_receipt': 'source_gate.json',
                      'comparison_receipt': 'full_content_QC_comparison.json'}
    for key, name in retained_names.items():
        if Path(member[key]) != folder / 'receipts' / name:
            raise RuntimeError('EXACT_PRIVATE_RETAINED_RECEIPT_ROUTE_REQUIRED')
    candidate_path = regular(spool / 'global_harmonization_candidate_receipt.json')
    column_path = regular(spool / 'columns/column_preparation_manifest.json')
    frozen_generated = {str(candidate_path): candidate_hash, str(column_path): columns_hash}
    unchanged(frozen_generated)
    candidate = json.loads(candidate_path.read_text())
    columns = json.loads(column_path.read_text())
    comparison_path = regular(member['comparison_receipt'])
    source_path = regular(member['source_gate_receipt'])
    worker_path = regular(root / 'receipts_v4' / (trait + '__compare.worker.json'))
    if sha(worker_path) != comparison_worker_hash:
        raise RuntimeError('FROZEN_COMPLETED_COMPARISON_WORKER_CHANGED')
    worker = json.loads(worker_path.read_text())
    expected_command = [plan['harmonization_python'], '-B', str(P / 'scripts/85_validate_core_large35_bounded_replay_v3.py'),
                        '--plan', str(plan_path), '--plan-sha', plan_hash,
                        '--trait', trait, '--mode', 'compare', '--out', str(comparison_path)]
    if (worker.get('status') != 'WORKER_COMPLETE_VERIFIED' or worker.get('plan_sha256') != plan_hash
        or worker.get('command') != expected_command or worker.get('returncode') != 0
        or worker.get('stop_reason') is not None or worker.get('metadata_errors')
        or worker.get('plan_unchanged') is not True or worker.get('owned_cleanup_verified') is not True
        or worker['process_group_teardown'].get('remaining_group_members')
        or worker['process_group_teardown'].get('cleanup_error')
        or worker['output_sha256'].get(str(comparison_path)) != sha(comparison_path)):
        raise RuntimeError('EXACT_COMPLETED_COMPARISON_WORKER_REQUIRED')
    for path in [worker_path, comparison_path, source_path]:
        if Path(str(path) + '.failure.json').exists() or Path(str(path) + '.failure.json').is_symlink():
            raise RuntimeError('FAILURE_ADDENDUM_PRESERVES_ALL_SCRATCH')
    comparison = json.loads(comparison_path.read_text())
    source = json.loads(source_path.read_text())
    raw = member['original_design']['raw']
    if (comparison.get('status') != 'EXACT_CORE_HARMONIZED_MUNGED_CONTENT_AND_QC_REPLAY_PASS'
        or source.get('status') != 'EXACT_CORE_RAW_SOURCE_GATE_PASS'
        or comparison.get('trait') != trait or source.get('trait') != trait
        or comparison.get('plan_sha256') != plan_hash or source.get('plan_sha256') != plan_hash
        or comparison.get('raw_sha256') != raw['sealed_verified_sha256']
        or source.get('raw_sha256') != raw['sealed_verified_sha256']
        or comparison.get('raw_bytes') != raw['bytes'] or source.get('raw_bytes') != raw['bytes']):
        raise RuntimeError('EXACT_FULL_CONTENT_QC_AND_SOURCE_PROOF_REQUIRED')
    h = comparison['harmonized_comparison']
    m = comparison['munged_comparison']
    rows = member['original_design']['harmonized_rows']
    template_rows = member['original_design']['historical_munge_log']['output_total_rows']
    if (h.get('status') != 'EXACT_ORDERED_DECOMPRESSED_HARMONIZED_MATCH'
        or h.get('rows') != [rows, rows] or h.get('expected_rows') != rows or h.get('unequal_rows') != 0
        or h.get('both_gzip_CRC_and_EOF_verified') is not True
        or len(h.get('decompressed_sha256', [])) != 2
        or any(not isinstance(d,str) or not re.fullmatch(r'[0-9a-f]{64}',d) for d in h['decompressed_sha256'])
        or h['decompressed_sha256'][0] != h['decompressed_sha256'][1]
        or m.get('status') != 'EXACT_FULL_DECOMPRESSED_TEMPLATE_MATCH'
        or m.get('rows_new') != template_rows or m.get('rows_archived') != template_rows
        or m.get('expected_template_rows') != template_rows or m.get('unequal_rows') != 0
        or m.get('both_gzip_CRC_and_EOF_verified') is not True
        or not isinstance(m.get('decompressed_sha256_new'),str) or not re.fullmatch(r'[0-9a-f]{64}',m['decompressed_sha256_new'])
        or not isinstance(m.get('decompressed_sha256_archived'),str) or not re.fullmatch(r'[0-9a-f]{64}',m['decompressed_sha256_archived'])
        or m['decompressed_sha256_new'] != m['decompressed_sha256_archived']
        or any(not isinstance(m.get(k),dict) or set(m[k])!={'finite_N_Z','missing_or_nonfinite_N_Z'}
            or any(type(v) is not int or v<0 for v in m[k].values()) or sum(m[k].values())!=template_rows
            for k in ['counts_new','counts_archived'])
        or m['counts_new'] != m['counts_archived']):
        raise RuntimeError('NESTED_EXACT_FULL_STREAM_AND_CRC_PROOFS_REQUIRED')
    expected_qc = {k: True for k in member['original_design']['historical_QC']['metadata']
                   if k not in {'infile', 'variant_map', 'liftover_chain'}}
    expected_qc['ordered_filter_steps'] = True
    if not expected_qc or comparison.get('QC_comparisons') != expected_qc:
        raise RuntimeError('ALL_EXACT_ORIGINAL_QC_AND_ORDERED_STEPS_REQUIRED')
    output_paths = [member[k] for k in ['harmonized', 'munged', 'harmonization_qc', 'source_gate_receipt']]
    output_hashes = comparison.get('output_identity_sha256', {})
    if set(output_hashes) != set(output_paths):
        raise RuntimeError('EXACT_FOUR_COMPLETED_COMPARISON_OUTPUT_BINDINGS_REQUIRED')
    unchanged(output_hashes)
    if h['compressed_sha256'][0] != output_hashes[member['harmonized']] or m['compressed_sha256_new'] != output_hashes[member['munged']]:
        raise RuntimeError('COMPARED_STREAM_AND_CURRENT_OUTPUT_IDENTITIES_DIFFER')
    if (candidate.get('status') != 'CANDIDATE_GLOBAL_HARMONIZATION_PRODUCED_NOT_SCIENTIFICALLY_ADMITTED'
        or candidate.get('original_code_sha256') != plan['dependencies_sha256'].get(member['literal_original_harmonize_command'][2])
        or candidate.get('source_sha256') != raw['sealed_verified_sha256']
        or candidate.get('expected_rows') != member['original_design']['harmonizer_input_rows']
        or candidate.get('rows_out') != rows or candidate.get('output_sha256') != output_hashes[member['harmonized']]
        or candidate.get('column_manifest_sha256') != columns_hash
        or columns.get('source_sha256_before') != raw['sealed_verified_sha256']
        or columns.get('source_sha256_after') != raw['sealed_verified_sha256']
        or columns.get('expected_rows') != member['original_design']['harmonizer_input_rows']):
        raise RuntimeError('FROZEN_GENERATED_METADATA_AND_FINAL_OUTPUT_BINDINGS_DIFFER')
    registered = {candidate_path, column_path, spool / 'first_surviving_SNP.sqlite3'}
    if '--variant-map' in member['literal_original_harmonize_command']:
        registered.add(spool / 'required_map_keys.sqlite3')
    pieces = [item for mapping in [columns['typed_pieces'], columns['preinference_token_pieces'],
                                  columns.get('coordinate_token_pieces', {})]
              for items in mapping.values() for item in items]
    pieces += candidate['post_liftover_frames'] + candidate['post_ordinary_QC_frames']
    for item in pieces:
        path = regular(item['path'])
        if (not path.is_relative_to(spool) or path in registered or sha(path) != item['sha256']
            or path.parent not in {spool / 'columns', spool / 'post_liftover', spool / 'post_ordinary_QC'}):
            raise RuntimeError('EXACT_UNIQUE_REGISTERED_GENERATED_PIECE_REQUIRED')
        registered.add(path)
    directories = {spool / name for name in ['columns', 'post_liftover', 'post_ordinary_QC']}
    observed_directories = set()
    files = set()
    for path in spool.rglob('*'):
        if path.is_symlink():
            raise RuntimeError('SPOOL_SYMLINK_PRESERVED_NO_DISPOSAL')
        if path.is_dir():
            observed_directories.add(path)
        elif path.is_file():
            files.add(path)
        else:
            raise RuntimeError('UNREGISTERED_SPECIAL_OBJECT_PRESERVED')
    if observed_directories != directories or not registered.issubset(files):
        raise RuntimeError('EXACT_REGISTERED_FILES_AND_DIRECTORIES_REQUIRED')
    sidecars = files - registered
    for path in sidecars:
        if not path.name.startswith('._') or path.with_name(path.name[2:]) not in registered | directories:
            raise RuntimeError('UNREGISTERED_FILE_OR_ORPHAN_SIDECAR_PRESERVED')
    inventory = {str(regular(path)): {'bytes': path.stat().st_size, 'sha256': sha(path)} for path in sorted(files)}
    proof = {**output_hashes, str(comparison_path): sha(comparison_path)}
    frozen_worker = {str(worker_path): comparison_worker_hash, **worker['output_sha256']}
    unchanged(frozen_worker)
    unchanged(frozen_generated)
    unchanged(proof)
    if sha(plan_path) != plan_hash:
        raise RuntimeError('CORE_PLAN_CHANGED_BEFORE_DISPOSAL')
    check_bindings(plan)
    copy_exact(candidate_path, member['bounded_receipt'], candidate_hash)
    copy_exact(column_path, member['column_manifest'], columns_hash)
    frozen = {'schema': 'owned_core_ephemeral_spool_inventory_v3', 'plan_sha256': plan_hash,
              'trait': trait, 'source_namespace': str(spool), 'scientific_success_proof_sha256': proof,
              'completed_comparison_worker_sha256': comparison_worker_hash,
              'registered_generated_paths': sorted(str(p) for p in registered),
              'registered_directories': sorted(str(p) for p in directories),
              'paired_ExFAT_sidecars': sorted(str(p) for p in sidecars), 'file_inventory': inventory,
              'generated_bytes': sum(v['bytes'] for v in inventory.values()), 'deletion_policy': POLICY,
              'copied_candidate_receipt_sha256': candidate_hash, 'copied_column_manifest_sha256': columns_hash}
    intended_inventory_hash=hashlib.sha256((json.dumps(frozen,indent=2)+'\n').encode()).hexdigest()
    write_new(member['spool_inventory'], frozen)
    if sha(regular(member['spool_inventory']))!=intended_inventory_hash:raise RuntimeError('PERSISTED_INVENTORY_DIFFERS_FROM_FROZEN_INTENT')
    retained_hashes = {member['bounded_receipt']: candidate_hash, member['column_manifest']: columns_hash,
                       member['spool_inventory']: intended_inventory_hash}
    # File contents and all directory entries are durable before any unlink.
    sync_directory(folder / 'receipts')
    unchanged(retained_hashes)
    unchanged(frozen_generated)
    unchanged(frozen_worker)
    unchanged(proof)
    if sha(plan_path) != plan_hash:
        raise RuntimeError('CORE_PLAN_CHANGED_AFTER_INTENT_PERSISTENCE')
    check_bindings(plan)
    deleted = 0
    for path, item in inventory.items():
        file = regular(path)
        if file.stat().st_size != item['bytes'] or sha(file) != item['sha256']:
            raise RuntimeError('GENERATED_SCRATCH_CHANGED_DURING_DISPOSAL')
        file.unlink()
        deleted += 1
    for directory in sorted(directories, reverse=True):
        if directory.is_symlink():
            raise RuntimeError('REGISTERED_DIRECTORY_CHANGED_DURING_DISPOSAL')
        directory.rmdir()
    spool.rmdir()
    sync_directory(spool.parent)
    unchanged(retained_hashes)
    unchanged(frozen_worker)
    unchanged(proof)
    if sha(plan_path) != plan_hash:
        raise RuntimeError('CORE_PLAN_CHANGED_AFTER_DISPOSAL')
    check_bindings(plan)
    record = {'schema': 'owned_core_ephemeral_spool_cleanup_v3',
              'status': 'OWN_GENERATED_SCRATCH_DISPOSED_AFTER_EXACT_FULL_CONTENT_QC_PASS',
              'completed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'plan_sha256': plan_hash, 'trait': trait, 'deleted_files': deleted,
              'deleted_bytes': frozen['generated_bytes'], 'spool_absent': not spool.exists(),
              'inventory_sha256': retained_hashes[member['spool_inventory']],
              'scientific_success_proof_sha256': proof,
              'completed_comparison_worker_sha256': comparison_worker_hash,
              'copied_candidate_receipt_sha256': candidate_hash,
              'copied_column_manifest_sha256': columns_hash, 'original_input_or_output_files_deleted': 0}
    cleanup_hash=hashlib.sha256((json.dumps(record,indent=2)+'\n').encode()).hexdigest()
    write_new(member['spool_cleanup_receipt'], record)
    if sha(regular(member['spool_cleanup_receipt']))!=cleanup_hash:raise RuntimeError('PERSISTED_CLEANUP_DIFFERS_FROM_FROZEN_INTENT')
    sync_directory(folder / 'receipts')
    unchanged(retained_hashes)
    unchanged({member['spool_cleanup_receipt']: cleanup_hash})
    unchanged(frozen_worker)
    unchanged(proof)
    if spool.exists() or spool.is_symlink() or sha(plan_path) != plan_hash:
        raise RuntimeError('FINAL_PRIVATE_DISPOSAL_STATE_CHANGED')
    check_bindings(plan)
    print(json.dumps({'trait': trait, 'status': record['status'],
                      'generated_bytes_disposed': record['deleted_bytes']}), flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--plan-sha', required=True)
    p.add_argument('--trait', required=True)
    p.add_argument('--expected-candidate-sha256', required=True)
    p.add_argument('--expected-columns-sha256', required=True)
    p.add_argument('--expected-comparison-worker-sha256', required=True)
    a = p.parse_args()
    run(a.plan, a.plan_sha, a.trait, a.expected_candidate_sha256,
        a.expected_columns_sha256, a.expected_comparison_worker_sha256)


if __name__ == '__main__':
    main()

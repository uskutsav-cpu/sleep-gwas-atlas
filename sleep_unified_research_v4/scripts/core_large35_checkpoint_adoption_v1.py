"""Narrow actual snoring QC adjudication and existing-result adoption checks."""
import json
from pathlib import Path

from extension_replay_common_v3 import sha

OLD_LABEL = 'N_eff below 50% of total (381,973.8) [CDG3]'
NEW_LABEL = 'sample size below 50% of configured effective N (381,973.8) [CDG3]'
RAW_SHA = '82eeb068648959afbbccf477745820a2118f88abdff8ff0076ebb47df21f2417'
HARMONIZED_SHA = '307c7d32a8d87fadf1d0c697b99c3cdd977351331cb0cbf4d74ad49de09d400b'
MUNGED_SHA = '44d46777b6011684c40e829056431c0a6a243cfeb433b6210d16f5050317b3d0'
FIRST8 = ['insomnia', 'sleepdur', 'shortsleep', 'longsleep', 'chronotype', 'sleepiness', 'napping', 'snoring']


def regular(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file() or any(p.is_symlink() for p in path.parents):
        raise RuntimeError('ADOPTED_EVIDENCE_MUST_REMAIN_REGULAR: ' + str(path))
    return path


def qc(path):
    metadata = {}; steps = []; inside = False
    for line in regular(path).read_text().splitlines():
        fields = line.split('\t')
        if fields == ['step', 'dropped', 'remaining']:
            inside = True
        elif inside and len(fields) == 3:
            steps.append(dict(reason=fields[0], dropped=int(fields[1]), remaining=int(fields[2])))
        elif len(fields) == 2:
            metadata[fields[0]] = fields[1]
    return metadata, steps


def semantic_snoring(trait, old, metadata, steps, columns, comparison):
    """One exact historical/current label pair, only in this actual N branch."""
    if trait != 'snoring' or columns['columns'].keys() & {'N', 'N_EFF', 'N_EFF_HALF', 'NCASE', 'NCONTROL'}:
        raise RuntimeError('ONLY_ACTUAL_SNORING_CONFIGURED_CONSTANT_N_BRANCH_ADMITTED')
    if metadata.get('sample_size_schema_mapping') != 'CONSTANT_N_EFF_FROM_MANIFEST':
        raise RuntimeError('EXACT_CONSTANT_N_BRANCH_PROOF_REQUIRED')
    if len(steps) != 15 or len(old['ordered_steps']) != 15:
        raise RuntimeError('EXACT15_SNORING_QC_STEPS_REQUIRED')
    expected = [dict(x) for x in old['ordered_steps']]
    if expected[13] != dict(reason=OLD_LABEL, dropped=0, remaining=7168629):
        raise RuntimeError('EXACT_ARCHIVED_SNORING_STEP14_REQUIRED')
    expected[13]['reason'] = NEW_LABEL
    if steps != expected:
        raise RuntimeError('SNORING_HAS_ADDITIONAL_QC_DIFFERENCES')
    if steps[14]['reason'] != 'sample-size mode: constant N_eff derived from config ncase/ncontrol':
        raise RuntimeError('EXACT_SNORING_SAMPLE_MODE_REQUIRED')
    omitted = {'infile', 'variant_map', 'liftover_chain'}
    checks = {k: metadata.get(k) == v for k, v in old['metadata'].items() if k not in omitted}
    if not all(checks.values()):
        raise RuntimeError('SNORING_SCIENTIFIC_METADATA_DIFFERS')
    expected_checks = {**checks, 'ordered_filter_steps': False}
    if comparison.get('status') != 'FAILED_PRESERVED_NO_AUTOMATIC_RETRY' or comparison.get('QC_comparisons') != expected_checks:
        raise RuntimeError('EXACT_PRESERVED_LABEL_ONLY_COMPARISON_FAILURE_REQUIRED')
    if comparison.get('trait') != trait or comparison.get('raw_sha256') != RAW_SHA or comparison.get('raw_bytes') != 258829716:
        raise RuntimeError('EXACT_SNORING_RAW_IDENTITY_PROOF_REQUIRED')
    h = comparison['harmonized_comparison']; m = comparison['munged_comparison']
    if (h.get('status') != 'EXACT_ORDERED_DECOMPRESSED_HARMONIZED_MATCH'
        or h.get('rows') != [7168629, 7168629] or h.get('expected_rows') != 7168629
        or h.get('unequal_rows') != 0 or h.get('decompressed_sha256') != [HARMONIZED_SHA] * 2
        or h.get('both_gzip_CRC_and_EOF_verified') is not True
        or m.get('status') != 'EXACT_FULL_DECOMPRESSED_TEMPLATE_MATCH'
        or [m.get(k) for k in ['rows_new', 'rows_archived', 'expected_template_rows']] != [1217311] * 3
        or m.get('unequal_rows') != 0 or m.get('decompressed_sha256_new') != MUNGED_SHA
        or m.get('decompressed_sha256_archived') != MUNGED_SHA
        or m.get('both_gzip_CRC_and_EOF_verified') is not True
        or m.get('counts_new') != dict(finite_N_Z=1179075, missing_or_nonfinite_N_Z=38236)
        or m.get('counts_archived') != m['counts_new']):
        raise RuntimeError('EXACT_ALREADY_COMPARED_SNORING_CONTENT_REQUIRED')
    return dict(status='SNORING_CONSTANT_N_SEMANTIC_QC_EQUIVALENCE_MATCHED_CONTENT_ADOPTED',
                alias_scope='SNORING_ONLY_CONFIGURED_CONSTANT_N_NO_SOURCE_N_COLUMNS',
                historical_label=OLD_LABEL, current_label=NEW_LABEL, step1=14,
                dropped=0, remaining=7168629, ncase=152302, ncontrol=256015,
                effective_N=4.0/(1.0/152302.0+1.0/256015.0), minimum_fraction=0.5,
                label_byte_identity=False, numerical_or_filter_threshold_change=False,
                matched_content_recalculated=False, old_failed_comparison_reclassified=False)


def adoption_gate(plan, *, current_content_hashes):
    """Authenticate small checkpoints now; current retained bytes at execution."""
    for path, expected in plan['checkpoint_metadata_sha256'].items():
        if sha(regular(path)) != expected:
            raise RuntimeError('PRESERVED_CORE_CHECKPOINT_METADATA_CHANGED: ' + path)
    old = json.loads(regular(plan['predecessor_plan']).read_text())
    failed = json.loads(regular(plan['predecessor_execution_receipt']).read_text())
    if (failed.get('status') != 'FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
        or failed.get('owned_cleanup_verified') is not True
        or failed.get('plan_sha256') != plan['predecessor_plan_sha256']
        or [x['trait'] for x in failed['completed_members']] != FIRST8[:7]
        or [x['trait_id'] for x in plan['adopted_members']] != FIRST8
        or [x['trait_id'] for x in old['members'][:8]] != FIRST8
        or len(plan['members']) != 27 or plan['member_count'] != 27
        or [x['trait_id'] for x in plan['members']] != [x['trait_id'] for x in old['members'][8:]]):
        raise RuntimeError('EXACT7_COMPLETED_PLUS_SNORING_PLUS27_CONTINUATION_REQUIRED')
    if plan['guard'] != old['guard'] or plan['deadline_anchor_utc'] != old['prepared_utc']:
        raise RuntimeError('UNCHANGED_CAPS_AND_CONSERVATIVE_ORIGINAL_DEADLINE_REQUIRED')
    def relocate(value):
        if isinstance(value, dict):return {k: relocate(v) for k, v in value.items()}
        if isinstance(value, list):return [relocate(v) for v in value]
        if isinstance(value, str) and value.startswith(old['private_namespace'] + '/'):
            return plan['private_namespace'] + value[len(old['private_namespace']):]
        return value
    if plan['members'] != [relocate(m) for m in old['members'][8:]] or plan['original35_order'] != [m['trait_id'] for m in old['members']]:
        raise RuntimeError('EXACT27_NAMESPACE_ONLY_COMMAND_TRANSFORMATION_REQUIRED')
    for item, prior in zip(plan['adopted_members'][:7], failed['completed_members']):
        if item['output_sha256'] != prior['output_sha256'] or item['member'] != old['members'][FIRST8.index(item['trait_id'])]:
            raise RuntimeError('EXACT_COMPLETED_RESULT_MAP_REQUIRED')
        comparison = json.loads(regular(item['member']['comparison_receipt']).read_text())
        if comparison.get('status') != 'EXACT_CORE_HARMONIZED_MUNGED_CONTENT_AND_QC_REPLAY_PASS' or comparison.get('plan_sha256') != plan['predecessor_plan_sha256']:
            raise RuntimeError('ADOPTED_COMPLETED_MEMBER_NOT_VERIFIED')
    snoring = plan['adopted_members'][7]; member = snoring['member']
    if member != old['members'][7]:raise RuntimeError('EXACT_ORIGINAL_SNORING_MEMBER_REQUIRED')
    metadata, steps = qc(member['harmonization_qc'])
    columns = json.loads(regular(snoring['columns_path']).read_text())
    comparison = json.loads(regular(member['comparison_receipt']).read_text())
    expected_semantics = semantic_snoring('snoring', member['original_design']['historical_QC'], metadata, steps, columns, comparison)
    addendum = json.loads(regular(plan['snoring_adjudication']).read_text())
    if addendum['semantic_decision'] != expected_semantics:
        raise RuntimeError('EXACT_SCOPED_SNORING_ADJUDICATION_REQUIRED')
    candidate = json.loads(regular(snoring['candidate_path']).read_text())
    if (candidate['column_manifest_sha256'] != snoring['output_sha256'][snoring['columns_path']]
        or candidate['output_sha256'] != snoring['output_sha256'][member['harmonized']]
        or snoring['output_sha256'][member['harmonized']] != comparison['harmonized_comparison']['compressed_sha256'][0]
        or snoring['output_sha256'][member['munged']] != comparison['munged_comparison']['compressed_sha256_new']):
        raise RuntimeError('SNORING_RETAINED_CONTENT_IDENTITY_DIFFERS_FROM_ACTUAL_COMPARISON')
    for path, digest in plan['adopted_worker_receipt_sha256'].items():
        record = json.loads(regular(path).read_text())
        if (sha(path) != digest or record.get('status') != 'WORKER_COMPLETE_VERIFIED'
            or record.get('plan_sha256') != plan['predecessor_plan_sha256']
            or record.get('command') != failed['expected_worker_commands'][path]
            or record.get('returncode') != 0 or record.get('stop_reason') is not None
            or record.get('metadata_errors') or record.get('owned_cleanup_verified') is not True
            or record['process_group_teardown'].get('remaining_group_members')
            or record['process_group_teardown'].get('cleanup_error')):
            raise RuntimeError('EXACT38_COMPLETED_CHECKPOINT_WORKERS_REQUIRED')
        if current_content_hashes:
            for output, expected in record['output_sha256'].items():
                if sha(regular(output)) != expected:raise RuntimeError('ADOPTED_WORKER_OUTPUT_CHANGED')
    if len(plan['adopted_worker_receipt_sha256']) != 38:
        raise RuntimeError('EXACT38_SUCCESSFUL_WORKERS_REQUIRED')
    if current_content_hashes:
        for item in plan['adopted_members']:
            for path, digest in item['output_sha256'].items():
                if sha(regular(path)) != digest:raise RuntimeError('ADOPTED_CONTENT_OR_METADATA_CHANGED')
    return dict(adopted_members=8, previously_completed=7, semantic_snoring=1,
                preserved_successful_workers=38, failed_predecessor_reclassified=False,
                current_retained_content_hashes_verified=current_content_hashes)

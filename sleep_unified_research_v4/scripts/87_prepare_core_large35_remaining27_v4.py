#!/usr/bin/env python3
"""Freeze the 27 unfinished members of the preserved original-core replay."""
import argparse
import json
import os
import shutil
from pathlib import Path

from extension_replay_common_v3 import sha, physical_mount, check_bindings

P = Path(__file__).resolve().parents[1]
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
V3 = SSD / 'core_pipeline/large35_bounded_replay_v3'
OUT = SSD / 'core_pipeline/large35_bounded_replay_v4_remaining27'
ADJUDICATION = P / 'reviews/independent_core35_snoring_QC_label_stop_review_v3_20261010.md'


def regular_hash(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file() or any(parent.is_symlink() for parent in path.parents):
        raise RuntimeError('EXACT_PRESERVED_REGULAR_EVIDENCE_REQUIRED: ' + str(path))
    return sha(path)


def qc_steps(path):
    rows = []
    in_steps = False
    for line in Path(path).read_text().splitlines():
        fields = line.split('\t')
        if fields == ['step', 'dropped', 'remaining']:
            in_steps = True
        elif in_steps and len(fields) == 3:
            rows.append({'reason': fields[0], 'dropped': int(fields[1]),
                         'remaining': int(fields[2])})
    return rows


def verify_v3(plan_path, receipt_path):
    plan_sha, receipt_sha = regular_hash(plan_path), regular_hash(receipt_path)
    plan, receipt = json.loads(plan_path.read_text()), json.loads(receipt_path.read_text())
    expected_order = [m['trait_id'] for m in plan['members']]
    completed = [m['trait'] for m in receipt.get('completed_members', [])]
    if (plan.get('member_count') != 35 or len(expected_order) != 35
        or receipt.get('status') != 'FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
        or receipt.get('estimator_calls') != 0
        or completed != expected_order[:7]
        or expected_order[7] != 'snoring'):
        raise RuntimeError('PRESERVED_V3_PREFIX_OR_STOP_STATE_DIFFERS')

    evidence = {str(plan_path): plan_sha, str(receipt_path): receipt_sha,
                str(ADJUDICATION): regular_hash(ADJUDICATION)}
    # Prove the seven already-completed traits remain intact without rerunning them.
    for item in receipt['completed_members']:
        for path, digest in item['output_sha256'].items():
            if regular_hash(path) != digest:
                raise RuntimeError('COMPLETED_V3_CORE_OUTPUT_CHANGED: ' + path)
            evidence[path] = digest
    for trait in completed:
        prefix = f"{V3}/receipts_v4/{trait}__"
        for path, digest in receipt['worker_receipt_sha256'].items():
            if path.startswith(prefix):
                if regular_hash(path) != digest:
                    raise RuntimeError('COMPLETED_V3_WORKER_RECEIPT_CHANGED: ' + path)
                evidence[path] = digest
    # The v3 snoring attempt stopped only because one QC label differed. Bind
    # the failed comparison itself and prove every other comparison is exact.
    snoring = next(m for m in plan['members'] if m['trait_id'] == 'snoring')
    compare_path = Path(snoring['comparison_receipt'])
    compare_sha = regular_hash(compare_path)
    comparison = json.loads(compare_path.read_text())
    old_steps = snoring['original_design']['historical_QC']['ordered_steps']
    new_steps = qc_steps(snoring['harmonization_qc'])
    old_reason = 'N_eff below 50% of total (381,973.8) [CDG3]'
    new_reason = 'sample size below 50% of configured effective N (381,973.8) [CDG3]'
    if (comparison.get('trait') != 'snoring'
        or comparison.get('plan_sha256') != plan_sha
        or comparison.get('status') != 'FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
        or comparison.get('error') != 'RuntimeError: COMPLETE_CORE_SERIALIZED_CONTENT_OR_QC_DIFFERS'
        or comparison.get('QC_comparisons', {}).get('ordered_filter_steps') is not False
        or any(v is not True for k, v in comparison.get('QC_comparisons', {}).items()
               if k != 'ordered_filter_steps')
        or comparison.get('harmonized_comparison', {}).get('status') != 'EXACT_ORDERED_DECOMPRESSED_HARMONIZED_MATCH'
        or comparison.get('harmonized_comparison', {}).get('rows') != [7168629, 7168629]
        or comparison.get('harmonized_comparison', {}).get('unequal_rows') != 0
        or comparison.get('harmonized_comparison', {}).get('both_gzip_CRC_and_EOF_verified') is not True
        or len(comparison.get('harmonized_comparison', {}).get('decompressed_sha256', [])) != 2
        or comparison['harmonized_comparison']['decompressed_sha256'][0] != comparison['harmonized_comparison']['decompressed_sha256'][1]
        or comparison.get('munged_comparison', {}).get('status') != 'EXACT_FULL_DECOMPRESSED_TEMPLATE_MATCH'
        or comparison.get('munged_comparison', {}).get('rows_new') != 1217311
        or comparison.get('munged_comparison', {}).get('rows_archived') != 1217311
        or comparison.get('munged_comparison', {}).get('unequal_rows') != 0
        or comparison.get('munged_comparison', {}).get('both_gzip_CRC_and_EOF_verified') is not True
        or comparison['munged_comparison'].get('decompressed_sha256_new') != comparison['munged_comparison'].get('decompressed_sha256_archived')
        or len(old_steps) != 15 or len(new_steps) != 15
        or old_steps[13] != {'reason': old_reason, 'dropped': 0, 'remaining': 7168629}
        or new_steps[13] != {'reason': new_reason, 'dropped': 0, 'remaining': 7168629}
        or any(old_steps[i] != new_steps[i] for i in range(15) if i != 13)):
        raise RuntimeError('SNORING_SINGLE_STEP14_LABEL_ADJUDICATION_DIFFERS')
    h, m = comparison['harmonized_comparison'], comparison['munged_comparison']
    expected_outputs = {snoring['harmonized']: h['compressed_sha256'][0],
                        snoring['munged']: m['compressed_sha256_new']}
    source_gate = json.loads(Path(snoring['source_gate_receipt']).read_text())
    if (source_gate.get('status') != 'EXACT_CORE_RAW_SOURCE_GATE_PASS'
        or source_gate.get('trait') != 'snoring'
        or source_gate.get('plan_sha256') != plan_sha
        or source_gate.get('raw_sha256') != comparison.get('raw_sha256')
        or source_gate.get('raw_bytes') != comparison.get('raw_bytes')):
        raise RuntimeError('SNORING_SOURCE_GATE_DIFFERS_FROM_FAILED_COMPARISON')
    expected_outputs[snoring['source_gate_receipt']] = regular_hash(snoring['source_gate_receipt'])
    expected_outputs[snoring['harmonization_qc']] = regular_hash(snoring['harmonization_qc'])
    for path, digest in expected_outputs.items():
        if regular_hash(path) != digest:
            raise RuntimeError('SNORING_COMPARED_OUTPUT_CHANGED: ' + path)
        evidence[path] = digest
    evidence[str(compare_path)] = compare_sha

    # Bind every successful snoring worker, its output/log/journal files, and
    # the failed compare worker that is deliberately absent from the v3 master.
    prefix = f"{V3}/receipts_v4/snoring__"
    for path, digest in receipt['worker_receipt_sha256'].items():
        if path.startswith(prefix):
            if regular_hash(path) != digest:
                raise RuntimeError('PRESERVED_SNORING_WORKER_RECEIPT_CHANGED: ' + path)
            evidence[path] = digest
            worker = json.loads(Path(path).read_text())
            if worker.get('plan_sha256') != plan_sha or worker.get('status') != 'WORKER_COMPLETE_VERIFIED':
                raise RuntimeError('PRESERVED_SNORING_SUCCESSFUL_WORKER_DIFFERS: ' + path)
            for output, output_sha in worker.get('output_sha256', {}).items():
                if regular_hash(output) != output_sha:
                    raise RuntimeError('PRESERVED_SNORING_WORKER_OUTPUT_CHANGED: ' + output)
                evidence[output] = output_sha
    failed_compare = V3 / 'receipts_v4/snoring__compare.worker.json'
    failed_worker_sha = regular_hash(failed_compare)
    failed_worker = json.loads(failed_compare.read_text())
    expected_command = receipt['expected_worker_commands'].get(str(failed_compare))
    if (failed_worker.get('status') != 'WORKER_FAILED_PRESERVED'
        or failed_worker.get('plan_sha256') != plan_sha
        or failed_worker.get('command') != expected_command
        or failed_worker.get('returncode') != 1
        or failed_worker.get('stop_reason') is not None
        or failed_worker.get('plan_unchanged') is not True
        or failed_worker.get('owned_cleanup_verified') is not True
        or failed_worker.get('process_group_teardown', {}).get('remaining_group_members')
        or failed_worker.get('output_sha256', {}).get(str(compare_path)) != compare_sha):
        raise RuntimeError('PRESERVED_SNORING_FAILED_COMPARATOR_RECEIPT_DIFFERS')
    evidence[str(failed_compare)] = failed_worker_sha
    for output, output_sha in failed_worker.get('output_sha256', {}).items():
        if regular_hash(output) != output_sha:
            raise RuntimeError('PRESERVED_SNORING_FAILED_WORKER_OUTPUT_CHANGED: ' + output)
        evidence[output] = output_sha
    adjudication = ADJUDICATION.read_text()
    if ('381,973.8' not in adjudication or '7,168,629' not in adjudication
        or '0 dropped' not in adjudication):
        raise RuntimeError('EXACT_SNORING_LABEL_ADJUDICATION_NOT_FOUND')
    return plan, receipt, expected_order, evidence


def replace_route(value, old, new):
    if isinstance(value, str):
        if value == old:
            return new
        if value.startswith(old + '/'):
            return new + value[len(old):]
        return value
    if isinstance(value, list):
        return [replace_route(v, old, new) for v in value]
    if isinstance(value, dict):
        return {k: replace_route(v, old, new) for k, v in value.items()}
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--v3-plan', type=Path, default=V3 / 'core_large35_bounded_replay_plan_v3.json')
    parser.add_argument('--v3-receipt', type=Path, default=V3 / 'core_large35_execution_receipt_v3.json')
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise RuntimeError('PRIOR_REVIEW_PLAN_PRESERVED_NO_OVERWRITE')
    # The reviewer copy lives in the writable project workspace. The real
    # frozen plan is emitted under OUT and requires the verified SSD mount.
    target_ssd = args.output.resolve().is_relative_to(OUT)
    if target_ssd:
        physical_mount()
    v3, receipt, order, excluded = verify_v3(args.v3_plan, args.v3_receipt)
    if OUT.exists() or OUT.is_symlink():
        raise RuntimeError('SUCCESSOR_NAMESPACE_ALREADY_EXISTS_NO_AUTOMATIC_RETRY')
    if target_ssd:
        internal_free = shutil.disk_usage('/System/Volumes/Data').free
        ssd_free = shutil.disk_usage(SSD).free
        if internal_free < v3['guard']['internal_floor_bytes'] or ssd_free < v3['guard']['SSD_floor_bytes']:
            raise RuntimeError('UNCHANGED_RESOURCE_GUARD_FAILS')
    members = []
    for source in v3['members'][8:]:
        member = json.loads(json.dumps(source))
        # Rebind only successor-owned paths and command arguments. Frozen historical
        # design records remain untouched as provenance for the original analysis.
        for key in ['prefilter_command', 'harmonize_command', 'munge_command']:
            member[key] = replace_route(member[key], str(V3), str(OUT))
        for key in ['ephemeral_spool', 'bounded_receipt', 'column_manifest', 'spool_inventory',
                    'spool_cleanup_receipt', 'source_gate_receipt', 'comparison_receipt',
                    'harmonized', 'harmonization_qc', 'munged', 'munged_prefix']:
            member[key] = replace_route(member[key], str(V3), str(OUT))
        members.append(member)
    if len(members) != 27 or [m['trait_id'] for m in members] != order[8:]:
        raise RuntimeError('EXACT_REMAINING27_ORDER_REQUIRED')
    dependencies = dict(v3['dependencies_sha256'])
    for path in [Path(__file__), P / 'scripts/88_run_core_large35_remaining27_v4.py',
                 P / 'scripts/85_validate_core_large35_bounded_replay_v3.py',
                 P / 'scripts/86_cleanup_core_bounded_spool_v3.py', ADJUDICATION]:
        dependencies[str(path)] = sha(path)
    result = dict(v3)
    result.update(
        schema='original_core_remaining27_bounded_raw_replay_plan_v4',
        source_v3_plan=str(args.v3_plan), source_v3_plan_sha256=sha(args.v3_plan),
        source_v3_receipt=str(args.v3_receipt), source_v3_receipt_sha256=sha(args.v3_receipt),
        preserved_completed_v3_members=order[:7], preserved_snoring_stop=True,
        snoring_adjudication_path=str(ADJUDICATION),
        snoring_adjudication_sha256=sha(ADJUDICATION),
        excluded_v3_evidence_sha256=excluded,
        members=members, member_count=27, private_namespace=str(OUT),
        dependencies_sha256=dependencies, archived_input_sha256={},
        scientific_membership_or_threshold_changes=False, estimator_calls=0,
        automatic_retry=False,
        family_reconciliation_required={
            'v3_completed_members': order[:7],
            'snoring_additive_adjudication_sha256': sha(ADJUDICATION),
            'successor_required_members': order[8:],
            'family_complete_only_if_all_35_reconcile': True,
            'successor_alone_is_not_full35_completion': True},
        execution_preconditions=(
            'Exact separate root admission and independent review; immutable original v3 artifacts; '
            'all190 baseline gate; original reviewed bounded harmonizer and stock-munge compatibility; '
            'one worker under shared mutex; unchanged resource guards; raw source identity and exact '
            'whole-output/QC comparisons; no estimator calls.'),
        predecessor_v3_receipt_status=receipt['status'])
    check_bindings(result)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if target_ssd:
        for member in members:
            for folder in ['prefilter', 'harmonized', 'munged', 'receipts']:
                (OUT / member['trait_id'] / folder).mkdir(parents=True, exist_ok=True)
        for folder in ['receipts_v4', 'logs_v4', 'tmp', 'cache']:
            (OUT / folder).mkdir(parents=True, exist_ok=True)
        if (shutil.disk_usage('/System/Volumes/Data').free < v3['guard']['internal_floor_bytes']
            or shutil.disk_usage(SSD).free < v3['guard']['SSD_floor_bytes']):
            raise RuntimeError('UNCHANGED_RESOURCE_GUARD_FAILS_AFTER_NAMESPACE_PREPARATION')
    with args.output.open('x') as f:
        f.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
        f.flush()
        os.fsync(f.fileno())
    print(json.dumps({'plan': str(args.output), 'sha256': sha(args.output),
                      'member_count': 27, 'members': [m['trait_id'] for m in members],
                      'workers_launched': 0, 'estimator_calls': 0,
                      'v3_completed_preserved': order[:7], 'snoring_rerun': False}, indent=2))


if __name__ == '__main__':
    main()

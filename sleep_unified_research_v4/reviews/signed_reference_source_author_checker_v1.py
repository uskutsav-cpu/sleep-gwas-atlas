#!/usr/bin/env python3
"""Metadata-only author closure; no source/runtime bodies, mutex or subprocess."""
import ast
import datetime
import hashlib
import json
import os
from pathlib import Path

P = Path(__file__).resolve().parents[1]
R = P / 'reviews'
PLAN = P / 'manifests/signed_reference_source_execution_plan_v1.json'
PLAN_SHA = '6d73f45ee3f73f7db7c0cc727d9b2d3e1af30d3ee27391cf1154f1f8ab8482b4'


def sha(p):
    p = Path(p)
    assert p.is_file() and not p.is_symlink(), ('nonregular', str(p))
    return hashlib.sha256(p.read_bytes()).hexdigest()


def consume(p, expected=None):
    h = sha(p)
    if expected is not None:
        assert h == expected, ('identity', str(p), h, expected)
    data = Path(p).read_bytes()
    assert hashlib.sha256(data).hexdigest() == h
    value = json.loads(data)
    assert sha(p) == h, ('after_parse', str(p))
    return value, h


def write_new(p, data):
    b = (json.dumps(data, sort_keys=True, indent=2) + '\n').encode()
    intended = hashlib.sha256(b).hexdigest()
    with Path(p).open('xb') as out:
        out.write(b)
        out.flush()
        os.fsync(out.fileno())
    assert sha(p) == intended
    return intended


def main():
    plan, _ = consume(PLAN, PLAN_SHA)
    deferred = plan['preparation_deferred_dependency_sha256']
    assert len(deferred) == 26
    fresh = {p: h for p, h in plan['dependency_sha256'].items() if p not in deferred}
    assert len(fresh) == 28
    before = {p: sha(p) for p in fresh}
    assert before == fresh
    script_names = ['signed_reference_source_prepare_v1.py', 'signed_reference_source_common_v1.py',
                    'signed_reference_source_controller_v1.py', 'signed_reference_source_worker_v1.py',
                    'signed_reference_decode_v1.py', 'signed_reference_ldsc_decoderB_v1.py']
    code = {str(P / 'scripts' / n): fresh[str(P / 'scripts' / n)] for n in script_names}
    for p in code:
        ast.parse(Path(p).read_text(), filename=p)
    assert plan['preparation_only'] is True and plan['execution_admitted'] is False
    for k in ['S_simulations_run', 'C_calibration_credit', 'estimator_calls']:
        assert plan[k] == 0
    assert plan['biological41Cov_release_admitted'] is False
    assert plan['resource_contract']['QC_controls_metadata_and_all_failed_output_cap_bytes'] == 512 * 1024**2
    bseal_path = R / 'independent_decoderB_runtime_identity_review_v1_seal.json'
    bseal, bseal_sha = consume(bseal_path, '7d7f388d4380632682f3dde920b32f93e2b7e454d90b2a72dac1afde3c7a0a51')
    bbindings = dict(bseal['file_sha256'])
    bbindings.update(bseal['bound_actual_profile_producers_receipts_sha256'])
    assert len(bseal['file_sha256']) == 6
    for p, h in bbindings.items():
        assert sha(p) == h
    root_path = P / 'logs/decoderB_runtime_independent_review_root_consumed_v1.json'
    root, root_sha = consume(root_path, '10f53c18adc86bb78b743b43951460b51050e832a5fc449812517bc897c16441')
    design_path = R / 'signed_ld_503eur_source_design_seal_v4_1.json'
    design, design_sha = consume(design_path, 'e63162e196fecece21f6114fe22d2563ff0857f61efb5ec4972af9c7cef60ead')
    for p, h in design['review_artifact_sha256'].items():
        assert sha(p) == h
    toy = {'regular_file_sha256': {}, 'symlink_literal': {}, 'physical_directories': [],
           'qualification': 'ALL INVENTED FIXTURES; NO ACTUAL SOURCE OR INDIVIDUAL DATA'}
    for name in ['signed_reference_source_fixture_v1', 'signed_reference_source_fixture_v1_2',
                 'signed_reference_source_fixture_v1_3', 'signed_reference_source_metadata_controls_v1']:
        directory = R / name
        for base, dirs, files in os.walk(directory, followlinks=False):
            toy['physical_directories'].append(str(Path(base)))
            for n in list(dirs) + files:
                p = Path(base) / n
                if p.is_symlink():
                    toy['symlink_literal'][str(p)] = os.readlink(p)
                elif p.is_file():
                    toy['regular_file_sha256'][str(p)] = sha(p)
    toy['physical_directories'].sort()
    toy_sha = write_new(R / 'signed_reference_source_invented_fixture_inventory_v1.json', toy)
    controls = {}
    for path, count in [('signed_reference_source_fixture_v1_3/fixture_receipt.json', 25),
                        ('signed_reference_source_metadata_controls_v1/metadata_control_receipt.json', 17)]:
        q, h = consume(R / path)
        assert q['case_count'] == count and all(x['status'] == 'PASS' for x in q['cases'])
        controls[str(R / path)] = {'sha256': h, 'case_groups': count}
    d = R / 'signed_reference_source_fixture_v1_3/invented_private/format/qc'
    paths = [d / 'static_selection.jsonl', d / 'J_ordered.jsonl']
    sizes = [p.stat().st_size for p in paths]
    assert sizes == [3822, 7144]
    assert all(len(p.read_bytes().splitlines()) == 17 for p in paths)
    combined = sum(sizes)
    storage = {'schema': 'INVENTED_JSONL_STORAGE_EXPOSURE_NOT_ACTUAL_SOURCE_FAILURE',
               'fixture_file_sha256': {str(p): sha(p) for p in paths},
               'invented_eligible_rows': 17, 'static_bytes': sizes[0], 'J_bytes': sizes[1],
               'static_bytes_per_row': sizes[0] / 17, 'J_bytes_per_row': sizes[1] / 17,
               'combined_bytes_per_row': combined / 17,
               'strict_QC_component_cap_bytes': 512 * 1024**2,
               'maximum_rows_at_this_fixture_serialization_before_other_QC_files': (512 * 1024**2 * 17) // combined,
               'hypothetical_all_master_rows_eligible': 1290028,
               'hypothetical_combined_bytes_rounded_up': (1290028 * combined + 16) // 17,
               'actual_eligible_row_count_observed': False,
               'actual_component_storage_RSS_and_elapsed_success_certified': False,
               'additional_coordinate_DB_exclusions_controls_metadata_need_storage': True,
               'strict_cap_unchanged_no_J_dropping_or_quota_increase': True,
               'remaining_review': 'INDEPENDENT_WHOLE_ACTUAL_PLAN_REVIEW; CONSIDER PROSPECTIVE COMPACT EXACTLY_EQUIVALENT SUCCESSOR BEFORE ADMISSION'}
    storage_sha = write_new(R / 'signed_reference_source_storage_exposure_v1.json', storage)
    after = {p: sha(p) for p in fresh}
    assert after == before == fresh and sha(PLAN) == PLAN_SHA
    receipt = {'schema': 'SOURCE_ONLY_AUTHOR_METADATA_PREPARATION_CLOSURE',
               'prepared_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'verdict': 'AUTHOR_CONTROLS_PASS_WITH_STORAGE_EXPOSURE_NO_SOURCE_EXECUTION_ADMISSION',
               'plan_path': str(PLAN), 'plan_sha256': PLAN_SHA,
               'new_six_code_sha256': code, 'fresh_nonbody_dependency_sha256_before': before,
               'fresh_nonbody_dependency_sha256_after': after,
               'inherited_deferred_identity_sha256_not_read_or_currently_hashed': deferred,
               'deferred_identity_count': 26, 'invented_control_receipts': controls,
               'invented_control_case_groups_total': 42,
               'fixture_inventory_sha256': toy_sha, 'storage_exposure_sha256': storage_sha,
               'decoderB_independent_runtime_seal_sha256': bseal_sha,
               'decoderB_independent_review_and_actual_metadata_sha256': bbindings,
               'decoderB_root_consumption_receipt_sha256': root_sha,
               'source_design_seal_sha256': design_sha,
               'source_design_six_artifact_sha256': design['review_artifact_sha256'],
               'runtime_profile_pass_scope': 'CURRENT_IDENTITY_ONLY; SYSTEM_AMBIENT_DYNAMIC_HISTORY_DECODER_SOURCE_SEPARATE',
               'original_fixture_runtime_pending_field_preserved': 'TRUE AT FIXTURE CREATION; SUBSEQUENT INDEPENDENT CURRENT_RUNTIME_IDENTITY_PASS CROSS_BOUND HERE',
               'actual_source_archive_genotype_reference_body_reads': 0,
               'actual_network_or_production_worker_or_mutex_operations': 0,
               'S_simulations_run': 0, 'C_calibration_credit': 0, 'estimator_calls': 0,
               'root_source_execution_admission': False,
               'whole_source_controller_independent_review_required': True}
    h = write_new(R / 'signed_reference_source_author_metadata_receipt_v1.json', receipt)
    print(json.dumps({'status': receipt['verdict'], 'receipt_sha256': h, 'case_groups': 42,
                      'new_frozen_code_count': 6, 'fresh_nonbody': 28, 'deferred_unread': 26,
                      'invented_file_count': len(toy['regular_file_sha256']), 'source_operations': 0}))


if __name__ == '__main__':
    main()

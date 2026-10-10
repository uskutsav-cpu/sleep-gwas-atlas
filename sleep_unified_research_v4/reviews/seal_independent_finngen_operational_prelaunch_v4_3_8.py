"""Seal narrow corrected prelaunch controls; no source/reference-body reads."""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

P = Path(__file__).resolve().parents[1]
R = P / 'reviews'

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def read(name):
    return json.loads((R / name).read_text())

def write_new(path, value):
    text = value if isinstance(value, str) else json.dumps(value, indent=2, allow_nan=False) + '\n'
    with Path(path).open('x', encoding='utf-8') as stream:
        stream.write(text)

plan_path = P / 'manifests/finngen_preprocessing_plan_v4_3_8.json'
fixed = {
    str(plan_path): '78f0faea56024be8eca65421bf655c0f7c224afe6f029ed3ca0cea22491781b8',
    str(P / 'scripts/58_run_finngen_feasibility_stage_v3_8.py'): 'b7167d99d50dc1ef97334a3a572fcb129e9eeaca07e4040a17c8a89defbf69c1',
    str(P / 'scripts/59_prepare_finngen_feasibility_stage_v3_8.py'): 'bce49437281170da5b70933082bd311a426e6ec36520cc21d42a628aa468f394',
    str(R / 'independent_finngen_operational_prelaunch_seal_v4_3_7.json'): '8e02e674783178c4765148cac16c00b114554b18c781af51276955c384e15af2',
}
assert all(sha(p) == digest for p, digest in fixed.items())
plan = json.loads(plan_path.read_text())
controls = read('independent_finngen_controller_contradiction_controls_receipt_v4_3_8.json')
binding = read('independent_finngen_prelaunch_binding_receipt_v4_3_8.json')
code = read('independent_finngen_code_correction_receipt_v4_3_8.json')
prior = read('independent_finngen_operational_prelaunch_seal_v4_3_7.json')
assert len(controls['controls']) == 23 and controls['all23_expected_controls_pass']
assert controls['candidate_unchanged'] and controls['no_actual_worker_body_reference_network_mutex']
assert binding['all_checks_pass'] and binding['dependencies_count'] == 641
assert len(binding['fresh_metadata_code_runtime_sha256']) == 594
assert len(binding['reference_body_hashes_bound_without_fresh_read']) == 47
assert binding['all476_rejection_artifacts_frozen_and_current']
assert binding['exact_170_donor_dependencies_retained']
assert code['resource_meter_AST_unchanged'] and code['scientific_worker_SHA_unchanged']
assert code['all_science_argv_equal_except_plan_epoch']
negative = [c for c in controls['controls'] if not c['expected_success']]
assert len(negative) == 21 and all(not c['execute_returned_success'] and not c['consumer_accepts_stage'] and c['PENDING_remains'] for c in negative)
corrected_labels = ['metadata_errors_nonempty', 'plan_unchanged_false', 'teardown_cleanup_error', 'wrong_command', 'wrong_returncode', 'stop_reason_nonempty', 'worker_mutation_after_parse']
assert all(next(c for c in controls['controls'] if c['label'] == label)['passed'] for label in corrected_labels)
assert all(sha(path) == digest for path, digest in prior['artifact_sha256'].items())

md_path = R / 'independent_finngen_operational_prelaunch_v4_3_8.md'
json_path = R / 'independent_finngen_operational_prelaunch_v4_3_8.json'
seal_path = R / 'independent_finngen_operational_prelaunch_seal_v4_3_8.json'
report = dict(
    schema='independent_finngen_operational_prelaunch_review_v4_3_8',
    completed_utc=datetime.now(timezone.utc).isoformat(),
    verdict='QUALIFIED_PRELAUNCH_PASS_CORRECTED_SEVEN_WORKER_RECEIPT_WITNESSES',
    prelaunch_pass=True, execution_admission_granted=False,
    actual_prepared_plan=str(plan_path), actual_prepared_plan_sha256=fixed[str(plan_path)],
    fixed_reviewed_identity_sha256=fixed,
    admitted_review_scope='Exact result-free preprocessing plan/controller correction, dependency/command/resource closure and bounded metadata fixtures only.',
    source_optimization_proof_inherited_unchanged=dict(token_cases=514, original_row_controls=34,
        Python='3.9.23', scientific_worker_sha256=plan['dependencies_sha256'][str(P / 'scripts/56_preprocess_finngen_insomnia_feasibility_v3.py')],
        row_control_sha256=plan['dependencies_sha256'][str(P / 'scripts/57_verify_finngen_row_controls_v3.py')],
        timing_not_actual_source_throughput=True),
    measured_controls=controls['controls'], all23_controls_expected_results=True,
    corrected_prior_fault_labels=corrected_labels,
    all21_negative_controls_reject_with_PENDING_preserved=True,
    worker_failure_addendum_control_rejected=True,
    exact_intended_command_and_full_receipt_oracle_reviewed=True,
    frozen_worker_SHA_before_first_parse_and_before_after_later_reads=True,
    hashes_after_output_consumption_checked=True,
    code_AST_unchanged_except_worker_read_correction_and_epoch=code,
    actual_plan_closure=dict(dependencies=641, all170_donor_dependencies_retained=True,
        all476_rejection_artifacts_in_dependencies=True, fresh_small_metadata_code_runtime=594,
        receipt_bound_reference_bodies_not_reread=47, science_inputs_argv_N_and_guards_unchanged=True,
        old9_failure_artifacts_current=True, old_FAILED_PENDING_empty_owned_group_preserved=True,
        current_outputs_and_root_admission_absent_at_binding_check=binding['new_attempt_outputs_and_admission_absent'],
        original190_completion_addendum_sha256=binding['original190_completion_addendum_sha256']),
    source_identity_receipt_binding_only=binding['source_binding'],
    reference_body_hash_bindings_without_reread=binding['reference_body_hashes_bound_without_fresh_read'],
    inherited_archival_symlink_exceptions=binding['exact_inherited_historical_symlink_bindings'],
    guard=binding['guard'], global_reservation_bytes=binding['global_reservation_bytes'],
    exact_single_preprocessing_command=binding['exact_single_preprocessing_command'],
    source_or_replication_admission=False,
    actual_scientific_workers=0, estimator_calls=0, network_calls=0,
    real_source_reference_body_parses=0, actual_heavy_mutex_acquisitions=0,
    root_admission_requirements=[
        'Bind exact current prepared plan, executor/preparer, this review/report/seal, all required independent review files and all641 operational dependencies in a separate root admission.',
        'Recheck all190 completion/current receipts, source/reference/runtime identities, physical SSD mount, current resource floors/caps and absence of prior result/PENDING/failure/seal in the new3_8 namespace.',
        'Queue exclusively through the same inherited-FD shared heavy mutex; never overlap an owned heavy worker. Preserve deferred signals and fully verified owned cleanup before lock release.',
        'Launch only the single frozen preprocessing command. Preserve old3_2 failure and all3_7 rejection artifacts; no automatic retry, overwrite, h2/rg fit or broader source/scientific admission follows from this review.',
    ],
    qualifications=[
        'Actual execute/Terminal2 logic was exercised with worker, monitor, lock, physical-mount, baseline and disk/time operations replaced by metadata mocks. Real process lifecycle and mutex exclusivity are inherited qualified controls.',
        'Receipt oracle checks were validated on the preserved seven counterexamples plus a genuine-schema healthy mock; this is not an exhaustive arbitrary malformed-JSON/type or concurrent-adversary certification.',
        'The source/reference bodies were not read or independently rehashed. Their full hashes are frozen prior-receipt bindings; actual preprocessing must enforce complete before/after hashes, CRC/EOF and ordered-output diagnostics.',
        'Native Python and w_hm3 inherited symlinks are allowed only at exact frozen literal/resolved targets. Runtime binary SHA was checked; reference-body SHA remained receipt-bound.',
        'Terminal2 is qualified by its privately owned immutable marker, deferred-signal and POSIX unlink contract; uncontrolled marker loss/replacement is outside the sealed guarantee.',
        'The exact stable-dictionary optimization proof and34 row controls remain valid, but toy lookup timings do not establish actual-source throughput or completion within two hours. No old complete derivative exists for direct content equivalence.',
        'Missing per-variant INFO/N, assumed effective N185146.70377332723, clinical/overlap/ascertainment questions and human source decisions remain unresolved. Computational success cannot admit source validity, independent replication, liability/population h2 or biological-null claims.',
    ])
write_new(json_path, report)
write_new(md_path, f'''# Exact FinnGen preprocessing 3_8 prelaunch review

**Qualified prelaunch PASS for the exact corrected preprocessing plan; root execution admission remains separate.** Plan SHA `{fixed[str(plan_path)]}`; executor SHA `{fixed[str(P / 'scripts/58_run_finngen_feasibility_stage_v3_8.py')]}`; preparer SHA `{fixed[str(P / 'scripts/59_prepare_finngen_feasibility_stage_v3_8.py')]}`. No actual pipeline worker, source/reference body parse, estimator, network operation or heavy mutex acquisition occurred.

All seven preserved 3_7 negative witnesses now reject: nonempty metadata_errors, false plan_unchanged, teardown cleanup_error, wrong command, nonzero returncode, nonempty stop_reason and bytes changed after first worker parse. The additional worker failure-addendum witness rejects. All 23 bounded controller controls passed their expected results: healthy execution/consumption and the no-retry first completion succeed; the 21 negative cases fail both execution and Terminal2 consumption while retaining private PENDING. The prior result/output/journal/source/stage/seal mutation and late resource/deadline/signal/failed-supplemental-write controls remain successful. Fixtures and original accepted 3_7 witnesses are preserved.

The new [read_fixed_worker helper]({P}/scripts/58_run_finngen_feasibility_stage_v3_8.py:123) checks the frozen regular worker SHA before/after parse and after consuming its output hashes, exact substituted intended command, complete status/plan, returncode zero, null stop, empty metadata_errors, plan_unchanged true, empty reaped group, empty cleanup error and owned-cleanup true, with a failure-addendum veto. The initial snapshot is now taken before first parse; final identity callbacks consume that same frozen SHA. These controls address receipt consumption defensively; they do not claim the genuine producer emitted the old contradictory receipts.

Namespace-normalized AST comparisons leave every top-level helper unchanged except execute's corrected worker consumption, and add only read_fixed_worker. The resource-meter AST is identical. All scientific argv and plan fields are unchanged except private plan/output epoch and prepared resource observations. The preparer now restricts its CLI to preprocessing, removing the unreachable h2 branch. The scientific56v3 worker and57v3 row controls retain their exact hashes; the sealed514 stable-dictionary equivalence cases and identical34 row controls are inherited without another biological read or throughput claim.

The real plan has641 dependency bindings. It retains all170 donor dependencies unchanged and binds every one of the476 rejection artifacts plus the exact rejection seal. Fresh hashes matched594 small metadata/code/runtime artifacts. The remaining47 reference-body identities were checked as frozen receipt bindings and path metadata without body reads; this review does not recompute raw/reference biological hashes. The acquired source remains a regular809,346,932-byte file with the same frozen generation/MD5/SHA binding. Exact inherited native-Python/w_hm3 symlinks retain their qualified literal/resolved-target identities. All nine old failure hashes, FAILED/PENDING, reaped worker and empty partial output remain preserved. New3_8 result/master/worker/PENDING/seal/root-admission paths were absent at the binding check.

The [frozen preprocessing detail]({P}/FROZEN_FINNGEN_PREPROCESSING_DETAIL_v1.md) and [feasibility protocol]({P}/FROZEN_FINNGEN_INSOMNIA_FEASIBILITY_PROTOCOL_v1.md) retain original build/allele/MAF/MHC/statistic filters, source order/eligible duplicate rules, deterministic direct-munged serialization and diagnostics without favorable-P selection. Constant effective N185146.70377332723 remains explicitly assumed; per-variant INFO and N are unverified. This review admits one result-free preprocessing command and zero fits. It admits no clinical/source qualification, independent replication, h2/rg result, liability/population-h2 or biological conclusion.

Operational safeguards remain one worker/BLAS thread, the inherited shared heavy-mutex FD,3GiB internal/5GiB SSD floors,2GiB owned RSS and full Finn namespace, two-hour deadline,300GiB campaign reservation, all190 prerequisite, deferred signals and verified owned cleanup. Mock controls do not independently re-certify actual process teardown or mutex exclusion. Terminal2 retains its sealed private immutable marker/deferred-signal/POSIX unlink boundary; uncontrolled marker loss is not covered. The actual-schema healthy fixture and measured fault corrections are tested; arbitrary malformed JSON/types or adversarial concurrent path replacement are not exhaustively certified.

Before launch, root must separately bind this exact plan/code/review/seal and required reviews, recheck all190 and all641 dependencies/source/reference/runtime/current resources/physical mount, ensure new3_8 result/terminal/failure paths remain unused, and queue under the same shared exclusive mutex after the current heavy worker is gone. Preserve old3_2 failures and3_7 rejection artifacts; no automatic retry or overwrite is authorized. A later observed-h2 plan requires the genuine completed3_8 preprocessing producer and separate prelaunch admission.
''')

files = {Path(path) for path in prior['artifact_sha256']}
files.update(Path(path) for path in fixed)
files.update([Path(__file__).resolve(), md_path, json_path])
for child in R.iterdir():
    if child.name.startswith('independent_finngen_') and 'v4_3_8' in child.name:
        if child.is_dir():
            files.update(p for p in child.rglob('*') if p.is_file() and not p.is_symlink())
        elif child.is_file() and not child.is_symlink():
            files.add(child)
assert all(p.is_file() and not p.is_symlink() for p in files)
file_map = {str(path): sha(path) for path in sorted(files)}
write_new(seal_path, dict(schema='independent_finngen_operational_prelaunch_seal_v4_3_8',
    verdict=report['verdict'], prelaunch_pass=True, execution_admission_granted=False,
    plan_sha256=fixed[str(plan_path)], executor_sha256=fixed[str(P / 'scripts/58_run_finngen_feasibility_stage_v3_8.py')],
    preserved_rejection_seal_sha256=fixed[str(R / 'independent_finngen_operational_prelaunch_seal_v4_3_7.json')],
    artifact_count=len(file_map), artifact_sha256=file_map, all23_expected_controls_pass=True,
    qualified_private_terminal_contract=True, scientific_source_admitted=False))
assert all(sha(path) == digest for path, digest in file_map.items())
print(json.dumps(dict(verdict=report['verdict'], artifact_count=len(file_map),
    sealed_files={str(path): sha(path) for path in [md_path, json_path, seal_path]})))

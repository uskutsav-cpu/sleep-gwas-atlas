"""Seal the completed independent extension-v7 metadata review, never execute it.

Reads only inspected code/plan/review metadata and this reviewer's own invented
metadata witnesses. It does not read inherited GWAS/reference/runtime payloads,
invoke candidate preparers/controllers, or acquire a campaign mutex.
"""
import datetime
import hashlib
import json
import os
from pathlib import Path
import resource
import time

P = Path(__file__).resolve().parents[1]
R = P / 'reviews'
S = P / 'scripts'
PREFIX = 'independent_extension_pipeline_v7'
REPORT = R / (PREFIX + '_prelaunch_review_v1.json')
MARKDOWN = R / (PREFIX + '_prelaunch_review_v1.md')
SEAL = R / (PREFIX + '_prelaunch_review_seal_v1.json')
RESOURCE_PLAN = R / (PREFIX + '_review_resource_plan_v1.json')
started = time.monotonic()


def bound():
    if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > (128 << 20):
        raise RuntimeError('INDEPENDENT_SEAL_RSS_EXCEEDS128MiB')
    if time.monotonic() - started > 180:
        raise RuntimeError('INDEPENDENT_SEAL_DEADLINE_EXCEEDS180s')


def sha(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise RuntimeError('SEAL_NONREGULAR_INPUT:' + str(path))
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(65536), b''):
            h.update(chunk)
    return h.hexdigest()


def save_new(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())


metadata_path = R / (PREFIX + '_metadata_receipt_v1.json')
controls_path = R / (PREFIX + '_controller_receipt_v1_1.json')
failure_path = R / (PREFIX + '_fixture_storage_failure_v1.json')
metadata = json.loads(metadata_path.read_text())
controls = json.loads(controls_path.read_text())
failed_fixture = json.loads(failure_path.read_text())
assert metadata['all_checks_pass'] is True
assert controls['all_expected_results'] is True and controls['controls_count'] == 27
assert len(controls['controls']) == 27 and all(c['passed'] is True for c in controls['controls'])
assert sum(c['metadata_worker_calls'] for c in controls['controls']) == 5214
assert controls['actual_worker_estimator_curl_mutex_body_runtime_read_calls'] == 0
assert controls['old402_1729_or190_reruns'] == 0
assert controls['fixture_bytes'] <= (64 << 20) and controls['max_RSS_bytes'] <= (128 << 20)
assert controls['elapsed_seconds'] <= 180
assert failed_fixture['exit_code'] == 1
assert failed_fixture['fixture_bytes'] > failed_fixture['declared_fixture_limit_bytes'] == (64 << 20)

current_metadata = dict(metadata['fixed_current_metadata_sha256'])
current_metadata.update(metadata['all24_new_dependency_metadata_hashes_current'])
for old_seal in metadata['current_preserved_author_rejection_and_old_preflight_seals']:
    current_metadata.update(old_seal['metadata_artifact_sha256'])
current_metadata[str(S / '47_run_extension_pipeline_replay_v6.py')] = '35ae2a55ff3523442c28d097df018771cc505e52ca29c818fe6a67b37a53a701'
current_metadata[str(S / 'terminal_commit_common_v2.py')] = '9ef14eda37d07b6a3442181820f8909bedd7764aae01fbcaafc966791dd5a7dd'
assert all(sha(path) == digest for path, digest in current_metadata.items())
bound()

save_new(RESOURCE_PLAN, dict(
    schema='independent_extension7_metadata_review_resource_plan_v1',
    control_run_limit_seconds=180, control_run_RSS_limit_bytes=128 << 20,
    control_run_fixture_limit_bytes=64 << 20, seal_limit_seconds=180,
    seal_RSS_limit_bytes=128 << 20, process_parallelism=1,
    real_worker_estimator_network_mutex_source_reference_runtime_payload_operations=0,
    preserved_first_fixture_storage_failure_bytes=failed_fixture['fixture_bytes'],
    successful_additive_fixture_bytes=controls['fixture_bytes'],
    original_failed_fixture_not_removed_or_rewritten=True,
    successful_harness_changes='Only compact invented operational-plan projection; production candidates and original caps unchanged.',
    seal_reads='Own invented metadata witness files plus inspected code/plan/review metadata; no inherited biological fixture/body or native rerun.'))

prerequisites = [
    'Separately issued exact-plan root execution admission, binding this independent review and author/preserved rejection evidence.',
    'Fresh full current verification of all518 frozen dependencies, including100 archived processed inputs;494 inherited identities are deferred here.',
    'Genuine successful complete100 raw8 producer before initial root admission/launch, including exact producer plan/admission/executor/Terminal2 and immutable per-source receipts/body identities.',
    'Current unchanged independently qualified full190 native completion and all acquisition source/version/HTTP/body/ownership/teardown gates.',
    'Current original3GiB internal/5GiB SSD floors,2GiB aggregate RSS,32GiB private output and300GiB campaign cap/ledger,1worker/BLAS1/shared mutex/deferred-signal owned cleanup.',
    'Actual success requires100 completed members, exact400 current worker identities/commands and full ordered decompressed original-template content comparisons under Terminal2; any failure remains preserved, no automatic retry.',
]
report = dict(
    schema='independent_actual_extension7_prelaunch_review_v1',
    created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    verdict='QUALIFIED_PREFLIGHT_PASS', execution_admission=False,
    actual_scientific_replay_or_source_clearance_certified=False,
    material_production_blockers_measured_in_review=[],
    author_review_independence='Separate reviewer of root-authored actual v7; author12 helper cases do not substitute for these27 whole-controller controls.',
    exact_plan_sha256=controls['actual_plan_sha256'], exact_executor_sha256=controls['reviewed_controller_sha256'],
    exact_preparer_sha256='c08a02c14d8290f02491e3a159e23c4e30628058aff3d47aa049a71e584af852',
    current_metadata_sha256=current_metadata,
    metadata_receipt=dict(path=str(metadata_path), sha256=sha(metadata_path)),
    controls_receipt=dict(path=str(controls_path), sha256=sha(controls_path)),
    resource_plan=dict(path=str(RESOURCE_PLAN), sha256=sha(RESOURCE_PLAN)),
    source_members=100, frozen_order=list(range(1, 101)), owned_commands=400,
    command_categories=dict(source_validation=100, historical_harmonization=100, historical_munging=100, full_content_comparison=100),
    estimator_commands=0, declared_raw_body_bytes=metadata['raw_declared_bytes'],
    declared_archived_processed_bytes=metadata['archived_processed_declared_bytes'],
    original_template_rows_per_source=1217311,
    all100_complete_scientific_members_unchanged_after_only_private_output_namespace_normalization=True,
    inherited_dependency_count=494, new_current_dependency_metadata_count=24, total_frozen_dependencies=518,
    inherited_dependency_or_archived_input_payloads_freshly_verified_here=False,
    old_v6_history='Preserved qualified13-artifact PASS plus later targeted6-artifact REJECT; v7 measured strict-oracle corrections reviewed additively without erasing either.',
    unchanged_AST=dict(top_level=metadata['unchanged_top_level_AST'], nested=metadata['unchanged_nested_AST']),
    controls=controls['controls'], controls_passed=27, fabricated_metadata_worker_events=5214,
    real_worker_estimator_network_mutex_operations=0, biological_source_reference_runtime_payload_reads=0,
    original402_controls_1729_fixture_payloads_190_numeric_rechecks=0,
    control_resources=dict(elapsed_seconds=controls['elapsed_seconds'], max_RSS_bytes=controls['max_RSS_bytes'], fixture_bytes=controls['fixture_bytes']),
    first_harness_failure=dict(path=str(failure_path), sha256=sha(failure_path),
        reviewer_fixture_storage_bound_failure=True, production_defect=False, exit_code=1,
        preserved_bytes=failed_fixture['fixture_bytes'], limit_bytes=64 << 20,
        corrected_only_in_additive_v1_1_private_harness=True),
    final_oracle_control_qualification=controls['qualifications'],
    terminal_qualification='Private immutable owned PENDING, deferred signals, normal POSIX unlink; no guarantee for uncontrolled marker deletion or effect-then-error mocked unlink. Read-only consumers must require committed exact receipt/current result identities.',
    provisional_success_is_not_committed_success='Some failing fixtures preserve a primary receipt with provisional ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS, but PENDING remains and both execution and Terminal2 consumer reject, including failed supplemental writes.',
    scientific_reproduction_qualification='Unchanged original harmonization, previously qualified existing Python3 header text compatibility for historical munging, original direct arguments. Actual full ordered decompressed SNP/A1/A2/Z/N/literal-missingness content and original-template identity remain future execution obligations; compressed gzip equality is not required.',
    original_filter_count_concordance='UNESTABLISHED_WHEN_QC_TSV_UNAVAILABLE',
    biological_claims_admitted=[], biological_covariance_power_novelty_independence_claims_admitted=[],
    campaign_guard=dict(internal_floor_bytes=3 << 30, SSD_floor_bytes=5 << 30,
        aggregate_RSS_limit_bytes=2 << 30, private_output_cap_bytes=32 << 30,
        global_cap_bytes=300 << 30, worker_count=1, BLAS_threads=1,
        worker_seconds=7200, family_seconds=345600),
    preserved_reservation=dict(reserved_total_bytes=320204645117, ceiling_bytes=322122547200,
        unallocated_margin_bytes=1917902083, new_unlimited_capacity_assumed=False),
    runner_initial_gate_qualification=metadata['qualification'],
    root_admission_prerequisites=prerequisites,
    reviewed_actual_output_attempt_pending_seal_master_absent=metadata['current_output_attempt_pending_seal_master_absent'],
)
save_new(REPORT, report)

text = f"""# Independent actual extension v7 prelaunch review

Verdict: **QUALIFIED_PREFLIGHT_PASS**, with no material production blocker measured in this bounded review. This is a separate review of root-authored v7, not execution admission, verification of current biological inputs, or certification of an actual100-source replay. Root admission remains absent at the inspected preparation state.

The exact actual plan is `{controls['actual_plan_sha256']}`; executor47 is `{controls['reviewed_controller_sha256']}`; preparer46 is `c08a02c14d8290f02491e3a159e23c4e30628058aff3d47aa049a71e584af852`. The [metadata checker]({R / (PREFIX + '_metadata_checker_v1.py')}) and [metadata receipt]({metadata_path}) verify these current identities,100 unique members in original order, and complete member dictionaries identical to v6 after only private output namespace normalization. The declared100 raw bodies total227,610,388,647 bytes and archived processed inputs872,644,553 bytes; these are frozen declarations, not new body reads. Four commands per source yield400 owned commands, including200 validation commands, and zero estimators.

Both historical reviews remain preserved: the earlier [v6 qualified PASS seal]({R / 'genomicsem_extension_pipeline_v6_preflight_seal.json'}) (`ce945e2289da6ce73e75899e2c54ac5977e0ea6def639a5e48c31e20b5358831`,13 metadata artifacts), and the later [targeted v6 REJECT seal]({R / 'extension_pipeline_v6_worker_consumer_addendum_v1_seal.json'}) (`35254ae81812f01121255ab34d7ce967d24e402b8a1b973007454c81a28c0d7c`,6 artifacts). The six consumer contradictions measured in the latter supersede unconditional reliance on the older PASS; they do not show that genuine normally guarded workers produced those contradictions. The current seven-artifact [author seal]({R / 'root_extension_pipeline_v7_author_preparation_v1_seal.json'}) is also current and retained. The author's12 helper cases are additional evidence; this independent review exercised the whole actual controller route on newly invented metadata.

The corrected [executor]({S / '47_run_extension_pipeline_replay_v7.py'}:32) enumerates exact400 worker paths and substituted commands. Its strict reader at line49 requires the receipt's frozen SHA before parsing and again after current output consumption, exact plan/command, integer return0, explicit null stop_reason, empty metadata_errors, plan_unchanged is True, owned cleanup, no remaining group members or cleanup_error, nonempty regular current output map and failure-addendum veto. Immediate snapshots are frozen before consumption at line259; final verification requires the exact400-key map at line217. The master persistence at line110 binds serialized intended SHA before writing and rejects postwrite drift. The unchanged nested limits/regular-output/validation AST and inherited top-level acquisition/source/mutex gates are recorded in the metadata receipt.

The [additive controller checker]({R / (PREFIX + '_controller_checker_v1_1.py')}) and [27-case receipt]({controls_path}) passed. A healthy fixture ran400 fabricated worker events, completed100 members and exact400 identities, committed Terminal2 and was consumed successfully. All26 negative fixtures rejected execution and terminal consumption while retaining private PENDING.

| Bounded controls | Count | Result |
| --- | ---: | --- |
| Healthy exact400/100 commit and read-only consumer |1|PASS|
| metadata_errors, plan_unchanged, cleanup_error, wrong command, nonzero return, nonnull stop, each immediate and final |12|All rejected|
| Receipt drift after parse/during output consumption, failure addendum, receipt/output symlinks, empty outputs, bool returncode |7|All rejected|
| Master postwrite mutation, missing final identity, old output drift, postpersist resource floor/signal, failed supplemental writes |6|All rejected with PENDING|
| Owned cleanup before mutex release |1|Cleanup ordering passed under mock ownership/lock|

The six final-oracle controls deliberately inject a contradictory private worker snapshot and matching mutable in-memory map at the mocked acquisition-family boundary, after healthy immediate consumption, to isolate final oracle coverage. This is not evidence that a genuine producer or external actor can mutate the in-memory master. Source/acquisition/full190/physical resource probes and process/lock operations are mocked; the candidate strict reader, command/cardinality checks, master intent and real private Terminal2 persistence are exercised. No actual worker, process group, heavy mutex, source/reference/runtime payload reader, network transfer, estimator or original402/1729/190 recheck was invoked.

Some negative fixtures preserve a provisional primary `ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS` string after postwrite failures. It does not establish success: PENDING remains and both the controller and Terminal2 consumer reject, even when supplemental failure writes are forced to fail. The result inherits the sealed Terminal2 contract of a private immutable owned PENDING marker, deferred signals and normal POSIX unlink. It makes no uncontrolled-marker-deletion or mocked effect-then-error unlink guarantee. Future consumers must require the committed terminal and recheck exact current receipts and results.

The first independent harness printed all27 expected control results but exited1 on its own final storage assertion:72,619,128 bytes exceeded64MiB. The [failure receipt]({failure_path}), original checker and complete witnesses are preserved. This was a reviewer fixture-size error, not a production finding. Additive v1_1 compacted only invented fixture plans; the real candidates, same128MiB RSS/64MiB fixture/180s bounds and full actual-plan comparison were retained. Its completed measurements are111.097 seconds,66,437,120 bytes peak RSS and46,120,885 fixture bytes. Both witness epochs are sealed rather than overwritten. The seal hashes only these own invented metadata files and inspected code/review/plan metadata; old biological fixture payloads remain untouched.

The scientific command matrix remains the original [Pan-UKBB harmonizer]({P.parent / 'discovery_extension/scripts/10_harmonize_panukbb.py'}) with exact local source, panel and recovered HM3 reference, plus unchanged frozen historical direct munging arguments from [original11]({P.parent / 'discovery_extension/scripts/11_munge_extension.sh'}). The munging executable retains the separately documented [existing Python3 header text compatibility qualification]({P / 'source_provenance/historical_munge_header_text_compatibility_v1.json'}); it is not claimed to be byte-identical pinned stock Python source. The unchanged comparison validator will require the entire1,217,311-row decompressed original template, exact ordered SNP/A1/A2/Z/N fields and literal missingness, and gzip CRC/EOF for both outputs. Compressed serialization equality is not required. No original harmonization QC TSV was established in the inherited bounded search, so original filter-count concordance remains `UNESTABLISHED_WHEN_QC_TSV_UNAVAILABLE`. Actual output reproduction and scientific/source qualifications remain unresolved until real execution evidence is reviewed.

All494 donor dependency bindings remain unchanged but explicitly deferred;24 new code/review metadata identities were freshly checked, yielding518 frozen dependencies. Preparation did not freshly verify source/reference/runtime bodies or100 archived processed files. All13 old PASS,6 rejection and7 author metadata-artifact hashes were checked;1729 old fixture payloads and prior body controls were not reread. The future root admission must freshly verify all518 current dependencies and archived inputs, and establish genuine complete100 raw8 producer evidence. The runner retains checkpoint overlap and enforces full100 at final seal; its initial gate alone does not enforce the stricter frozen precondition requiring complete100 before root admission/launch. Root must enforce that ordering separately. The reported active90/100 acquisition state is not eligible for that admission.

Original guards remain: one heavy worker/BLAS1, shared mutex across each source's four commands and verified cleanup, unlocked checkpoint waiting, deferred SIGINT/TERM/HUP,3GiB internal/5GiB SSD floors,2GiB aggregate RSS,2h worker/96h family,32GiB pipeline cap and300GiB campaign meter. The preserved ledger reserves320,204,645,117 of322,122,547,200 bytes, leaving1,917,902,083 bytes; no extra unlimited capacity is assumed. Current full190 completion, source/version/HTTP/strong body hash/ownership gates and independently issued exact-plan root admission remain mandatory. Any stop preserves failure evidence and requires review; no automatic retry is admitted.

No new h2/rg outcome, biological association, covariance, power, novelty or independent replication claim follows from this review. This pass admits only the measured prospective controller correction under its stated private-owned evidence and inherited methods boundaries.
"""
with MARKDOWN.open('x') as stream:
    stream.write(text)
    stream.flush()
    os.fsync(stream.fileno())
bound()

# Hash only local invented witnesses; neither inherited494 payloads nor old1729
# fixture bodies are traversed. Do not follow the private symlink fault fixtures.
files = dict(current_metadata)
symlinks = {}
for root in [R / (PREFIX + '_control_fixtures_v1'), R / (PREFIX + '_control_fixtures_v1_1')]:
    assert root.is_dir() and not root.is_symlink()
    for path in root.rglob('*'):
        if path.is_symlink():
            symlinks[str(path)] = os.readlink(path)
        elif path.is_file():
            files[str(path)] = sha(path)
        bound()
for path in [metadata_path, controls_path, failure_path, REPORT, MARKDOWN, RESOURCE_PLAN,
    Path(__file__), R / (PREFIX + '_metadata_checker_v1.py'),
    R / (PREFIX + '_controller_checker_v1.py'), R / (PREFIX + '_controller_checker_v1_1.py')]:
    files[str(path)] = sha(path)
assert all(files[path] == digest for path, digest in controls['fixture_regular_sha256'].items())
assert all(symlinks[path] == target for path, target in controls['fixture_symlinks'].items())
bound()
seal = dict(schema='independent_actual_extension7_review_seal_v1',
    verdict='QUALIFIED_PREFLIGHT_PASS', execution_admission=False,
    file_sha256=files, regular_file_count=len(files), symlink_literal_targets=symlinks,
    exact_plan_sha256=controls['actual_plan_sha256'], exact_executor_sha256=controls['reviewed_controller_sha256'],
    metadata_and_invented_fixture_files_only=True, inherited494_payloads_and_old1729_fixtures_not_rehashed=True,
    all27_control_results_bound=True, original_failed_and_additive_successful_fixture_epochs_preserved=True,
    seal_max_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    seal_elapsed_seconds=time.monotonic() - started)
save_new(SEAL, seal)
for path, digest in files.items():
    assert sha(path) == digest
    bound()
assert all(Path(path).is_symlink() and os.readlink(path) == target for path, target in symlinks.items())
print(json.dumps(dict(verdict=seal['verdict'], execution_admission=False,
    report=str(REPORT), report_sha256=sha(REPORT), markdown=str(MARKDOWN), markdown_sha256=sha(MARKDOWN),
    seal=str(SEAL), seal_sha256=sha(SEAL), regular_files=len(files), symlinks=len(symlinks),
    elapsed_seconds=time.monotonic() - started, max_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)))

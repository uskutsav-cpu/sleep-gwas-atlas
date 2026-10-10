"""Seal conditional draft review, deliberately distinct from actual prelaunch."""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

P = Path(__file__).resolve().parents[1]
R = P / 'reviews'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write_new(path, value):
    text = value if isinstance(value, str) else json.dumps(value, indent=2, allow_nan=False) + '\n'
    with Path(path).open('x', encoding='utf-8') as stream:
        stream.write(text)

controls_path = R / 'independent_finngen_h2_candidate_controls_receipt_v4_3_9.json'
closure_path = R / 'independent_finngen_h2_candidate_closure_receipt_v4_3_9.json'
controls = json.loads(controls_path.read_text())
closure = json.loads(closure_path.read_text())
assert controls['all22_expected_results'] and len(controls['controls']) == 22
assert controls['reviewed_inputs_unchanged'] and closure['all_checks_pass']
assert all(sha(path) == digest for path, digest in closure['fixed_candidate_producer_consumer_sha256'].items())
assert all(closure['actual_h2_plan_admission_stage_absent'].values())
md_path = R / 'independent_finngen_h2_candidate_review_v4_3_9.md'
json_path = R / 'independent_finngen_h2_candidate_review_v4_3_9.json'
seal_path = R / 'independent_finngen_h2_candidate_review_seal_v4_3_9.json'
report = dict(
    schema='independent_finngen_h2_conditional_candidate_review_v4_3_9',
    completed_utc=datetime.now(timezone.utc).isoformat(),
    verdict='CONDITIONAL_CANDIDATE_PASS', conditional_candidate_pass=True,
    actual_plan_prepared=False, actual_plan_prelaunch_review_complete=False,
    h2_execution_admission_granted=False, scientific_source_admitted=False,
    fixed_draft_producer_consumer_sha256=closure['fixed_candidate_producer_consumer_sha256'],
    scope='Prospective code/binding correction for one observed-scale standalone FinnGen R13 F5_INSOMNIA h2 diagnostic; no rg, liability conversion or new scientific/source admission.',
    all22_pure_helper_controls_expected_results=True, controls=controls['controls'],
    code_and_producer_binding_checks=closure['AST_and_binding_checks'],
    prospective_prior_preprocessing=dict(plan_path=str(P / 'manifests/finngen_preprocessing_plan_v4_3_8.json'),
        plan_sha256='78f0faea56024be8eca65421bf655c0f7c224afe6f029ed3ca0cea22491781b8',
        stage_path=str(P / 'logs/finngen_preprocessing_stage_receipt_v4_3_8.json'),
        derivative_suffix='pipeline_replay_v3_8/derived/insomnia.sumstats.gz',
        result_suffix='pipeline_replay_v3_8/derived/insomnia.preprocessing.json',
        completed_stage_and_full_result_map_and_Terminal2_required=True,
        stage_result_identity_snapshots_before_parse_and_after=True,
        no_actual_producer_success_inferred_from_code_or_root_preprocessing_admission=True),
    qualified_controls_inherited=dict(
        exact3_8_worker_reader_and_execute_AST=True, prior_producer_gate3_6_AST=True,
        worker_snapshot_before_parse_and_later_frozen_reads=True,
        exact_intended_argv_returncode_stop_metadata_plan_cleanup_oracle=True,
        original23_control_execution_receipt_not_rerun=True,
        original514_source_semantic_cases_and34_row_controls_not_rerun=True,
        native190_numerical_reproduction_not_reaudited=True,
        final_consumer60v5_unchanged_sha256='e9e85f9a79d6cc07da4da9b41b3a0d6bd9b95a969b4165d31a0217e0950719e2'),
    preserved_review_history=closure['current_preserved_seals'],
    captured_data_verification_required_after_fit=[
        'Exactly one kind=h2 job/estimate/ordered intersection, original 200 blocks, pinned stock estimator and capture closure, matching actual argv and before/after input hashes.',
        'Native Terminal2 and current preprocessing Terminal2, exact capture/log/.delete/.part_delete/worker membership and hashes in completed result maps, regular paths and parents before/after consumption.',
        'Captured ordered SNP count and SNP/Z/N digests, one-annotation fields, finite full-precision h2/SE/intercept/lambda/mean-chi-square and current source qualifications.',
        'Stock total and partition delete vectors agree with capture; centered h2/intercept SE arithmetic and fixed relative1e-12/absolute1e-15 tolerances. Intercept uses captured native delete array because stock exports no separate intercept-delete file.',
        'Durable diagnostic JSON/TSV report pair and read-only require_verified_report veto for PENDING/failure/drift, including repeated current preprocessing gate after report-seal consumption.',
    ],
    actual_h2_plan_admission_stage_absent_at_review=closure['actual_h2_plan_admission_stage_absent'],
    actual_scientific_workers=0, fits=0, raw_reference_body_reads=0, network_calls=0, actual_heavy_mutex_calls=0,
    prospective_preparer_or_executor_calls=0, arithmetic_verifier_calls=0,
    qualifications=[
        'Tiny prior producer/capture identities are invented fixtures. Only candidate require_admission/result_gate and private Terminal2 metadata helpers were exercised; no biological estimate or arithmetic calibration was produced.',
        'The actual h2 plan and actual completed preprocessing producer do not yet exist as reviewed evidence. Candidate code cannot grant execution admission or substitute for a new exact producer-bound prelaunch review.',
        'This review rehashes prior small review/control artifacts and exact drafts; it does not freshly hash/decompress the biological raw source/reference bodies or redo the190 estimators.',
        'Stage completion gates enforce command/result scope and source bindings. Detailed stock-vector, final-intersection, scalar arithmetic and report-pair verification remains the unchanged postfit60v5 consumer duty; no such actual result exists yet.',
        'Terminal2 relies on its documented private immutable-owned PENDING, deferred-signal and POSIX unlink contract; arbitrary concurrent replacement/uncontrolled marker deletion is outside the guarantee.',
        'Effective N185146.70377332723 is assumed constant, with source INFO and per-variant N unverified. Computational/source h2 success cannot resolve clinical construct, chronicity/observation horizon, ascertainment, source independence or replication validity.',
        'Necessary diagnostic h2 Z>=4/intercept<=1.2 thresholds remain unchanged and cannot certify population h2 or scientific source validity; no liability prevalence or novel rg/BH/outcome test is admitted.',
        'One worker/BLAS1, inherited shared mutex, original3/5GiB floors,2GiB owned RSS/full Finn root, two-hour deadline,300GiB campaign meter and all190 prerequisites remain inherited safeguards, not a new live-process test.',
    ],
    remaining_root_actions=[
        'Wait for genuine admitted3_8 preprocessing execution to finish successfully with exact Terminal2 and complete result/hash/CRC/ordered-output proofs; retain any failure without automatic retry.',
        'Execute the new59 preparer only after that success, freezing actual producer metadata/terminal/input identities, all closure dependencies and the exact one observed-scale h2 command in the new namespace.',
        'Request a separate independent exact actual-plan prelaunch review and root operational admission. Candidate review filenames deliberately do not satisfy actual required prelaunch filenames.',
        'Queue exclusively under the shared heavy-worker mutex; after a fit, execute unchanged60v5 verification and require its verified durable report pair before consuming scientific diagnostics.',
    ])
write_new(json_path, report)
drafts = closure['fixed_candidate_producer_consumer_sha256']
write_new(md_path, f'''# Prospective FinnGen observed-h2 3_9 candidate review

**CONDITIONAL_CANDIDATE_PASS.** Exact 58 draft SHA `{drafts[str(P / 'scripts/58_run_finngen_feasibility_stage_v3_9.py')]}` and 59 draft SHA `{drafts[str(P / 'scripts/59_prepare_finngen_feasibility_stage_v3_9.py')]}` passed this narrow code/helper review. This is deliberately distinct from actual-plan prelaunch review and execution admission. No actual h2 plan, admission or stage receipt existed at the metadata check. No preparer/executor, scientific worker, fit, raw/reference body read, network operation, real heavy mutex or arithmetic verifier ran.

The active prospective producer is the exact actual 3_8 preprocessing plan SHA `78f0faea56024be8eca65421bf655c0f7c224afe6f029ed3ca0cea22491781b8`, its 3_8 stage/result/derivative and 59 producer. Old 3_2 FAILED/PENDING evidence and the unlaunched 3_6 donor remain preserved. The new preparer fixes the 3_8 plan SHA and snapshots prior plan/stage/result hashes before parsing, rechecks after parsing and before its final freeze, requires exact completed preprocessing status/result map and Terminal2, and binds prior admission/executor/preparer/seal identities. Preprocessing root admission alone is not completed-producer evidence.

All 22 pure-helper cases passed: healthy admission/result consumption succeeded; the other 21 rejected. Cases covered prior PENDING/failure/FAILED status, dependency/review/plan/executor/guard/scope drift, two jobs/estimates/intersections, mismatching input hashes, rg/population/sample-prevalence options, capture bytes changed immediately after parsing, and prior derivative/metadata drift during capture consumption. The new result SHA is frozen before parse and checked after consumption; the repeated prior gate rejects the measured later prior-input drift. Private fixture captures contain invented values and were not interpreted as biological results.

The worker reader and execute AST are exactly the qualified 3_8 implementation; its full argv/return/stop/metadata/plan/cleanup oracle, fixed preparse worker SHA, failure-addendum veto and final immutable output proofs therefore inherit the sealed 23 controls. The prior preprocessing gate is unchanged from the qualified 3_6 donor. The h2 scientific job update, source identity/effective-N expressions and baseline/runtime route are unchanged. The only active source-path change points to completed 3_8, and the output epoch changes to observed_h2_v3_9. The CLI admits only observed_h2. No 514-case source/34-row or original 190-job reproduction was repeated.

The preparer preserves all 1,175 earlier rejection artifact entries, the 62-artifact final 60v5 consumer review and the 943-artifact accepted 3_8 review; their exact seals and current small-artifact hashes were independently checked. Metadata rehashing does not freshly certify the biological source/reference bodies. The actual source assumptions remain those in the [feasibility protocol]({P}/FROZEN_FINNGEN_INSOMNIA_FEASIBILITY_PROTOCOL_v1.md) and [preprocessing detail]({P}/FROZEN_FINNGEN_PREPROCESSING_DETAIL_v1.md): constant effective N = 185146.70377332723 is assumed, INFO/per-variant N are unverified, and clinical/overlap/ascertainment/source-admission decisions remain unresolved.

The [unchanged 60v5 consumer]({P}/scripts/60_verify_finngen_h2_diagnostic_v5.py:88), SHA `e9e85f9a79d6cc07da4da9b41b3a0d6bd9b95a969b4165d31a0217e0950719e2`, remains required after a genuine fit. It verifies the exact one observed-h2 command with 200 blocks and no rg/liability options, current native/prior Terminal2, capture/log/total-delete/partition-delete/worker closure, final ordered SNP count and SNP/Z/N digests, captured args, stock total/partition vectors and scalar SE arithmetic. Intercept arithmetic uses the captured native delete array; stock has no separate intercept-delete export. Fixed relative 1e-12 and absolute 1e-15 arithmetic tolerances remain unchanged. Detailed result certification belongs to that postfit consumer; this draft review supplies no native result or sampling calibration.

Its durable report-pair and [read-only consumer]({P}/scripts/60_verify_finngen_h2_diagnostic_v5.py:58) require absent private PENDING/failure evidence and current native/prior/report hashes, repeating the prior gate after report-seal consumption. Both report and native Terminal2 retain the private immutable-owned marker/deferred-signal/POSIX unlink boundary; uncontrolled marker deletion or arbitrary concurrent replacement is outside the guarantee. No source admission, population h2/liability conversion, independent replication, new rg/BH/outcome test or biological-null claim follows from stage/report success. Original h2 Z>=4 and intercept<=1.2 boundaries are necessary diagnostics only.

Root must wait for actual successful 3_8 preprocessing, freeze the actual completed producer-bound h2 plan with this prospective code, obtain a separate exact actual-plan prelaunch review and root operational admission, and queue under the existing exclusive shared mutex. One worker/BLAS1, 3 GiB internal/5 GiB SSD floors, 2 GiB owned RSS/full Finn namespace, two-hour deadline, 300 GiB campaign reservation, all190 prerequisite and deferred owned cleanup remain unchanged. Preserve failures and all older artifacts; no automatic retry is admitted.
''')
files = {Path(path) for path in closure['fixed_candidate_producer_consumer_sha256']}
files.update([Path(__file__).resolve(), md_path, json_path])
for row in closure['current_preserved_seals']:
    seal = Path(row['seal'])
    files.add(seal)
    files.update(Path(path) for path in json.loads(seal.read_text())['artifact_sha256'])
for path in R.iterdir():
    if path.name.startswith('independent_finngen_h2_candidate_') and 'v4_3_9' in path.name:
        if path.is_dir():
            files.update(child for child in path.rglob('*') if child.is_file() and not child.is_symlink())
        elif path.is_file() and not path.is_symlink():
            files.add(path)
assert all(path.is_file() and not path.is_symlink() for path in files)
artifact_map = {str(path): sha(path) for path in sorted(files)}
write_new(seal_path, dict(schema='independent_finngen_h2_conditional_candidate_seal_v4_3_9',
    verdict='CONDITIONAL_CANDIDATE_PASS', actual_plan_prelaunch_review_complete=False,
    execution_admission_granted=False, all22_pure_helper_controls_pass=True,
    artifact_count=len(artifact_map), artifact_sha256=artifact_map,
    fixed_draft_identity_sha256=drafts, qualified_private_terminal_contract=True))
assert all(sha(path) == digest for path, digest in artifact_map.items())
print(json.dumps(dict(verdict='CONDITIONAL_CANDIDATE_PASS', artifact_count=len(artifact_map),
    sealed_files={str(path): sha(path) for path in [md_path, json_path, seal_path]})))

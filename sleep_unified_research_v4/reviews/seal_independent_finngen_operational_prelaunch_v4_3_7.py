"""Seal the exact measured rejection; metadata and already-created tiny fixtures only."""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

P = Path(__file__).resolve().parents[1]
R = P / 'reviews'

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def write_new(path, value):
    text = value if isinstance(value, str) else json.dumps(value, indent=2, allow_nan=False) + '\n'
    with Path(path).open('x', encoding='utf-8') as handle:
        handle.write(text)

plan_path = P / 'manifests/finngen_preprocessing_plan_v4_3_7.json'
expected = {
    str(plan_path): 'd75fb3f96795f38080081be89da10be04c82d74aa67d11821a7eac0d6e7a8046',
    str(P / 'scripts/56_preprocess_finngen_insomnia_feasibility_v3.py'): '9ae9ec603d9b8558d30ec48c69e6658586853986c2c6e6a474b4fbdb20ce7302',
    str(P / 'scripts/57_verify_finngen_row_controls_v3.py'): '35a969111c8dcdfca4e7f5bea619f936cbe2eda32b263dc13f359b8e4d199e12',
    str(P / 'scripts/58_run_finngen_feasibility_stage_v3_7.py'): '94183216b0db8e8be4c1d89b6b882b0b7453fdf97e3677fc7130c15e5760ef9c',
    str(P / 'scripts/59_prepare_finngen_feasibility_stage_v3_7.py'): 'c8f1e6136257f56c00339e500fc0feeba9e90d323be8167855909e34a7f1b11a',
    str(P / 'logs/finngen_preprocessing_performance_stop_v4_3_2.json'): '77e331f828d43669742b7d4b03320e2630acebffb30506a77791a843e38a185b',
}
for name, digest in expected.items():
    assert sha(name) == digest, name

def receipt(name):
    return json.loads((R / name).read_text())

eq = receipt('independent_finngen_rsid_equivalence_receipt_v4_3_7.json')
rows = receipt('independent_finngen_original34_row_recheck_v4_3_7.json')
initial = receipt('independent_finngen_controller_controls_receipt_v4_3_7.json')
negative = receipt('independent_finngen_controller_contradiction_controls_receipt_v4_3_7.json')
binding = receipt('independent_finngen_prelaunch_binding_receipt_v4_3_7.json')
assert eq['all514_equality_checks_pass'] and eq['original34_row_controls_rerun_identical']
assert rows['all_pass'] and len(rows['checks']) == 34
assert initial['all15_expected_controls_pass'] and len(initial['controls']) == 15
assert binding['all_checks_pass'] and len(binding['fresh_metadata_code_runtime_sha256']) == 123
assert len(binding['reference_body_hashes_bound_without_fresh_read']) == 47
assert len(negative['controls']) == 8 and negative['controls'][0]['passed']
unexpected = [c for c in negative['controls'] if not c['passed']]
assert len(unexpected) == 7
assert all(c['execute_returned_success'] and c['consumer_accepts_stage'] and not c['PENDING_remains'] for c in unexpected)

md_path = R / 'independent_finngen_operational_prelaunch_v4_3_7.md'
json_path = R / 'independent_finngen_operational_prelaunch_v4_3_7.json'
seal_path = R / 'independent_finngen_operational_prelaunch_seal_v4_3_7.json'
report = {
    'schema': 'independent_finngen_operational_prelaunch_review_v4_3_7',
    'completed_utc': datetime.now(timezone.utc).isoformat(),
    'verdict': 'REJECTED_CURRENT_CONTROLLER_REQUIRES_ADDITIVE_CORRECTION',
    'prelaunch_pass': False,
    'execution_admission_granted': False,
    'actual_prepared_plan': str(plan_path),
    'actual_prepared_plan_sha256': expected[str(plan_path)],
    'candidate_and_old_failure_receipt_sha256': expected,
    'source_optimization_semantic_equivalence_pass': True,
    'exact_only_scientific_code_change': 'set(RSID.findall(value)).intersection(hm3) -> {rsid for rsid in RSID.findall(value) if rsid in hm3}',
    'stable_ordinary_dict_token_cases': 514,
    'original_row_controls_repeated_identically': 34,
    'native_python': eq['python'],
    'native_probe_peak_RSS_bytes': eq['max_RSS_bytes'],
    'native_probe_seconds': eq['elapsed_seconds'],
    'dictionary_timing_rows': eq['timings'],
    'initial_controller_controls_expected_results': initial['controls'],
    'additive_receipt_contradiction_controls': negative['controls'],
    'material_accepted_witnesses': [c['label'] for c in unexpected],
    'findings': [
        {'code': 'WORKER_RECEIPT_ORACLE_INCOMPLETE', 'file': str(P / 'scripts/58_run_finngen_feasibility_stage_v3_7.py'), 'lines': [156, 157, 171, 172],
         'observed': 'Six independently malformed complete receipts were committed: metadata_errors nonempty, plan_unchanged false, teardown cleanup_error truthy, wrong command, nonzero returncode, nonempty stop_reason.',
         'required_correction': 'Require exact substituted intended command; returncode == 0; stop_reason is None; plan_unchanged is True; metadata_errors == []; owned cleanup is True; remaining group == []; no cleanup error. Reject missing or malformed required fields. Apply oracle at first consumption and each current-identity callback.'},
        {'code': 'WORKER_HASH_SNAPSHOT_AFTER_PARSE', 'file': str(P / 'scripts/58_run_finngen_feasibility_stage_v3_7.py'), 'lines': [156, 171],
         'observed': 'Appending a byte immediately after initial worker JSON parse was accepted: the changed bytes, rather than the parsed snapshot, were subsequently frozen.',
         'required_correction': 'Require regular path/parents and freeze exact worker SHA before initial read/parse; compare after parse, before use, and against this frozen SHA before/after every later parse and terminal callback.'},
    ],
    'producer_versus_consumer_qualification': 'Genuine owned-worker producer normally rejects these contradictions. Toy malformed receipts demonstrate defensive controller-consumption gaps; no genuine 3_7 worker was run or shown to emit contradictory metadata.',
    'metadata_binding_checks_pass': True,
    'fresh_small_metadata_code_runtime_dependencies': 123,
    'reference_hash_bindings_without_body_reread': 47,
    'raw_source_binding_not_fresh_rehash': binding['source_binding'],
    'inherited_archival_symlink_identity_exceptions': binding['exact_inherited_historical_symlink_bindings'],
    'old_failed_attempt_hashes_current': binding['preserved9_failed_attempt_hashes_current'],
    'old_FAILED_PENDING_group_cleanup_preserved': binding['old_failed_PENDING_and_empty_owned_group_preserved'],
    'guard': binding['guard'],
    'global_reservation_bytes': binding['global_reservation_bytes'],
    'scope': {'preprocessing_commands': 1, 'h2_fits': 0, 'rg_fits': 0, 'actual_scientific_workers': 0, 'actual_heavy_mutex_acquisitions': 0, 'source_reference_body_parses': 0, 'network_calls': 0},
    'scientific_admission': False,
    'qualifications': [
        'Stable ordinary dictionary with string keys only; not an arbitrary mutating custom mapping.',
        'Toy dictionary timings establish elimination of reference-size scaling in membership lookup, not actual-source throughput or completion within the two-hour deadline.',
        'No complete old preprocessing output exists for full-content comparison: old failure retains a zero-byte derivative and no completed source CRC/EOF proof.',
        'Source SHA/MD5/generation/bytes and reference SHAs are frozen receipt bindings. This review does not freshly read or hash their biological bodies.',
        'Inherited runtime and w_hm3 symlinks are admitted only at exact frozen literal/resolved targets; fresh runtime binary SHA checked, w_hm3 body SHA receipt-bound.',
        'Controller tests exercise actual execute/Terminal2 logic with worker, monitor, lock, physical mount, baseline and disk/time replaced by metadata mocks; owned live-process cleanup and actual mutex exclusivity are inherited qualified controls.',
        'Terminal2 requires its documented privately owned immutable marker, deferred-signal and POSIX unlink contract; uncontrolled marker loss is not certified.',
        'Constant effective N and missing INFO/per-variant N remain assumptions. No source, clinical, independent-replication, population-h2, liability or biological inference is admitted.',
        'Root must freeze an additive corrected controller/plan, independently recheck these witnesses, then separately issue operational admission before any real worker. No automatic retry or overwrite is authorized.',
    ],
}
write_new(json_path, report)
write_new(md_path, f'''# Exact FinnGen preprocessing 3_7 prelaunch review

**REJECTED for execution pending an additive controller correction.** Exact prepared plan SHA `{expected[str(plan_path)]}` and controller SHA `{expected[str(P / 'scripts/58_run_finngen_feasibility_stage_v3_7.py')]}` were reviewed. No real pipeline worker, source/reference body parse, estimator, network operation or heavy mutex acquisition was performed. This review grants no execution admission.

The source optimization is semantically admissible for the actual stable ordinary HM3 dictionary. The new membership set equals the old intersection set, including duplicate-collapse and single/multiple/empty disposition. Only that expression changed in 56v3. All 514 independent token cases agreed under the declared Python 3.9.23, and the original 34 row controls were repeated identically. With a toy 400,000-key dictionary, empty-token median time changed from 0.0103690 s to 0.0000008542 s; six timing cases agreed. The probe used 89,702,400 bytes peak RSS. These are bounded membership timings, not source throughput or a guarantee of completion within two hours. The old failed attempt has no completed derivative for content equivalence certification.

The first 15 metadata-only controller controls passed their expected results, including result snapshot mutation, output/journal/worker drift after persistence, exact stage/seal payload mutation, late resource/deadline/signal failures, failed supplemental writes, source CRC/INFO contradiction and no automatic retry. The additive eight-case worker-receipt controls then exposed seven accepted negative witnesses:

| Witness | Intended result | Measured result |
|---|---|---|
| metadata_errors nonempty | reject | execute and Terminal2 consumer accepted |
| plan_unchanged false | reject | execute and Terminal2 consumer accepted |
| teardown cleanup_error truthy | reject | execute and Terminal2 consumer accepted |
| command differs from intended argv | reject | execute and Terminal2 consumer accepted |
| returncode nonzero | reject | execute and Terminal2 consumer accepted |
| stop_reason nonempty | reject | execute and Terminal2 consumer accepted |
| worker bytes changed immediately after first parse | reject | changed bytes frozen; execute and consumer accepted |

The healthy additive case passed. Every negative case produced a committed stage with PENDING removed; all fixtures are retained. At [58v3_7 line 156]({P}/scripts/58_run_finngen_feasibility_stage_v3_7.py:156), the current worker oracle checks only complete status, plan hash, empty remaining group and owned-cleanup truthiness. At [line 171]({P}/scripts/58_run_finngen_feasibility_stage_v3_7.py:171), it parses before freezing the worker SHA. The existing result receipt uses a stronger hash-before-parse gate; the worker receipt needs the same exact-byte binding.

An additive controller must require exact intended command after plan/lock-FD substitution, returncode zero, null stop_reason, plan_unchanged true, empty metadata_errors, owned cleanup true, empty group members and absent/empty cleanup error. Required fields must reject missing or malformed values. Freeze the regular worker file SHA before its first read/parse, recheck after parse and before use, and compare every later read with that frozen SHA before/after parsing. Apply the full oracle both immediately after the worker and at final current-identity gates. Genuine worker-producer error logic normally blocks such contradictions; these witnesses do not show a genuine worker produced malformed evidence.

The exact real-plan metadata audit passed: 170 dependency bindings comprise 123 freshly hashed small code/metadata/runtime files and 47 receipt-bound reference bodies whose contents were not read. The source is a regular 809,346,932-byte file; frozen acquisition generation, MD5 and SHA agree with the plan. Raw/reference body hashes were not freshly recomputed by this review. Two inherited archival symlinks are qualified only at their exact literal and resolved targets (native Python and w_hm3); the initial overly strict checker failure and its explicit corrected checker are preserved. All nine old failure-artifact hashes remain current, with FAILED/PENDING, reaped worker, empty derivative and no success seal preserved. Current 3_7 output/admission paths were absent at the metadata check. No retry was launched.

The [frozen preprocessing detail]({P}/FROZEN_FINNGEN_PREPROCESSING_DETAIL_v1.md) and [source-feasibility protocol]({P}/FROZEN_FINNGEN_INSOMNIA_FEASIBILITY_PROTOCOL_v1.md) retain the same allele/build/MAF/MHC/supplied-statistic filters, first eligible source-order duplicates, deterministic direct-munged schema and formatting, diagnostics without favorable-P filtering, and assumed constant effective N = 185146.70377332723. INFO and per-variant N remain unverified; minimum retention is a necessary diagnostic only. No new scientific threshold, historical scope, clinical/source clearance, independent replication, h2/rg result, liability conversion or biological conclusion is admitted.

One worker/BLAS thread, inherited shared mutex FD, 3 GiB internal/5 GiB SSD floors, 2 GiB aggregate RSS and full Finn namespace, two-hour deadline, 300 GiB campaign reservation, all-190 prerequisite, deferred signals and verified owned cleanup remain frozen. This bounded mock review does not re-certify real process teardown. Terminal2 remains qualified by its private immutable marker and deferred-signal/POSIX unlink contract. Root must preserve this rejection, prepare a new namespace/controller/plan, obtain exact independent correction review and separately issue operational admission before execution. No automatic retry is permitted.
''')

prefixes = [
    'independent_finngen_rsid_equivalence_probe_v4_3_7',
    'independent_finngen_rsid_equivalence_receipt_v4_3_7',
    'independent_finngen_original34_row_recheck_v4_3_7',
    'independent_finngen_controller_controls_v4_3_7',
    'independent_finngen_controller_controls_receipt_v4_3_7',
    'independent_finngen_controller_contradiction_controls_v4_3_7',
    'independent_finngen_controller_contradiction_controls_receipt_v4_3_7',
    'independent_finngen_prelaunch_binding_checker_v4_3_7',
    'independent_finngen_prelaunch_binding_checker_correction_v4_3_7',
    'independent_finngen_prelaunch_binding_receipt_v4_3_7',
]
files = {Path(name) for name in expected}
files.update([Path(__file__).resolve(), md_path, json_path, P / 'logs/finngen_row_controls_v4_3_7.json',
              P / 'FROZEN_FINNGEN_PREPROCESSING_DETAIL_v1.md', P / 'FROZEN_FINNGEN_INSOMNIA_FEASIBILITY_PROTOCOL_v1.md'])
for child in R.iterdir():
    if any(child.name.startswith(prefix) for prefix in prefixes):
        if child.is_dir():
            files.update(p for p in child.rglob('*') if p.is_file() and not p.is_symlink())
        elif child.is_file() and not child.is_symlink():
            files.add(child)
file_map = {str(p): sha(p) for p in sorted(files)}
write_new(seal_path, dict(schema='independent_finngen_prelaunch_rejection_seal_v4_3_7',
    verdict=report['verdict'], prelaunch_pass=False, execution_admission_granted=False,
    plan_sha256=expected[str(plan_path)], artifact_count=len(file_map), artifact_sha256=file_map,
    artifacts_regular=True, qualified_private_terminal_contract=True,
    preserved_negative_witnesses=report['material_accepted_witnesses']))
assert all(sha(name) == digest for name, digest in file_map.items())
print(json.dumps({'verdict': report['verdict'], 'artifact_count': len(file_map),
    'sealed_files': {str(p): sha(p) for p in [md_path, json_path, seal_path]}}))

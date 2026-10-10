"""Actual plan/code/review metadata; 494 prior body identities remain deferred."""
import ast
import hashlib
import json
from pathlib import Path

P = Path(__file__).resolve().parents[1]
R = P / 'reviews'
S = P / 'scripts'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OLD = SSD / 'extension_pipeline_replay_v6/extension_pipeline_replay_plan_v6.json'
NEW = SSD / 'extension_pipeline_replay_v7/extension_pipeline_replay_plan_v7.json'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

fixed = {
    str(OLD): 'd0fd223d468f26583bfd2bd477a51548646b57e91275d7733eb1ae814b94b798',
    str(NEW): 'ad0c08bc0b0a62c8a1bc070221689f2733b0eefe8f25de8cc1bd8689a81b1ba5',
    str(S / '46_prepare_extension_pipeline_replay_v7.py'): 'c08a02c14d8290f02491e3a159e23c4e30628058aff3d47aa049a71e584af852',
    str(S / '47_run_extension_pipeline_replay_v7.py'): '6abd57c9d743ea1c33b7ccbed0f88c51a7af0972db24655284f3e279fa483f00',
    str(R / 'extension_pipeline_v6_worker_consumer_addendum_v1_seal.json'): '35254ae81812f01121255ab34d7ce967d24e402b8a1b973007454c81a28c0d7c',
    str(R / 'genomicsem_extension_pipeline_v6_preflight_seal.json'): 'ce945e2289da6ce73e75899e2c54ac5977e0ea6def639a5e48c31e20b5358831',
    str(R / 'root_extension_pipeline_v7_author_preparation_v1_seal.json'): 'bc3a68353750b1be2549a0e514875192f0230693b278c20a29542327adf5a089',
}
assert all(sha(path) == digest for path, digest in fixed.items())
old, new = json.loads(OLD.read_text()), json.loads(NEW.read_text())
def norm(value):
    if isinstance(value, str): return value.replace(str(NEW.parent) + '/', str(OLD.parent) + '/')
    if isinstance(value, list): return [norm(x) for x in value]
    if isinstance(value, dict): return {k: norm(v) for k, v in value.items()}
    return value
assert old['members'] == norm(new['members'])
assert len(new['members']) == len({m['extension_trait_id'] for m in new['members']}) == 100
assert [m['index'] for m in new['members']] == list(range(1, 101))
assert new['total_owned_commands'] == 400 and new['estimator_calls'] == 0 and new['template_rows'] == 1217311
assert new['guard'] == old['guard']
assert new['scientific_membership_or_threshold_changes'] is False and new['automatic_retry'] is False
assert new['full_decompressed_munged_byte_and_field_identity_required'] is True
assert new['compressed_harmonized_or_munged_byte_identity_required'] is False
assert len(old['dependencies_sha256']) == 494
assert new['inherited_dependency_sha256_not_freshly_verified_at_preparation'] == old['dependencies_sha256']
assert all(new['dependencies_sha256'][k] == v for k, v in old['dependencies_sha256'].items())
assert len(new['dependencies_sha256']) == 518
assert new['fresh_source_reference_runtime_body_verification_at_preparation'] is False
assert new['actual_execution_admission_requires_full_current_dependencies_including_archived_inputs'] is True
science_keys = ['acquisition_plan', 'acquisition_plan_sha256', 'acquisition_execution_binding', 'acquisition_family_terminal_contract',
    'acquisition_operational_identity_by_origin', 'archived_input_sha256', 'baseline_dependency_sha256', 'baseline_execution_plan_sha256',
    'baseline_input_sha256', 'baseline_relocated_dependency_sha256', 'baseline_support_package', 'checkpoint_receipt_binding_policy',
    'environment', 'python', 'ldsc_dir', 'original_190_jobs', 'panel', 'reference', 'reference_sha256', 'reference_prefix', 'w_hm3',
    'historical_munge_compatibility_provenance', 'historical_munge_compatibility_provenance_sha256', 'historical_munge_qualification',
    'original_filter_count_concordance', 'original_harmonization_QC_bounded_search', 'reservation_arithmetic',
    'global_resource_ledger_path', 'global_resource_ledger_sha256']
assert all(old[k] == new[k] for k in science_keys)
fresh = {path: digest for path, digest in new['dependencies_sha256'].items() if path not in old['dependencies_sha256']}
assert len(fresh) == 24 and all(sha(path) == digest for path, digest in fresh.items())
history = []
for filename, count, key in [
    ('extension_pipeline_v6_worker_consumer_addendum_v1_seal.json', 6, 'file_sha256'),
    ('root_extension_pipeline_v7_author_preparation_v1_seal.json', 7, 'file_sha256'),
    ('genomicsem_extension_pipeline_v6_preflight_seal.json', 13, 'artifacts'),
]:
    path = R / filename
    receipt = json.loads(path.read_text())
    mapping = receipt[key]
    assert len(mapping) == count
    hashes = {k: v['sha256'] if isinstance(v, dict) else v for k, v in mapping.items()}
    assert all(sha(k) == v for k, v in hashes.items())
    history.append(dict(path=str(path), seal_sha256=sha(path), metadata_artifacts=count, metadata_artifact_sha256=hashes,
        old_fixture_payloads_not_rehashed=True))
oldtree, newtree = ast.parse((S / '47_run_extension_pipeline_replay_v6.py').read_text()), ast.parse((S / '47_run_extension_pipeline_replay_v7.py').read_text())
oldfunc = {n.name: n for n in oldtree.body if isinstance(n, ast.FunctionDef)}
newfunc = {n.name: n for n in newtree.body if isinstance(n, ast.FunctionDef)}
unchanged = [name for name in oldfunc if name not in ['run', 'seal_master', 'main']]
assert all(ast.dump(oldfunc[name], include_attributes=False) == ast.dump(newfunc[name], include_attributes=False) for name in unchanged)
def nested(fn, name):
    return next(n for n in ast.walk(fn) if isinstance(n, ast.FunctionDef) and n.name == name)
for name in ['limits', 'regular_output_hashes', 'validation']:
    assert ast.dump(nested(oldfunc['run'], name), include_attributes=False) == ast.dump(nested(newfunc['run'], name), include_attributes=False)
assert set(newfunc) - set(oldfunc) == {'read_fixed_worker', 'expected_worker_commands'}
absent_paths = [NEW.parent / ('pipeline_' + name + '_v7.json') for name in ['attempt', 'pending', 'terminal_seal', 'execution_receipt']]
absent = {str(path): not (path.exists() or path.is_symlink()) for path in absent_paths}
assert all(absent.values())
raw_bytes = sum(m['raw_bytes'] for m in new['members'])
archived_bytes = sum(m['archived_munged_bytes'] for m in new['members'])
source_ledger = [dict(index=m['index'], trait=m['extension_trait_id'], raw_sha256=m['raw_sha256'], raw_md5=m['raw_md5'],
    raw_bytes=m['raw_bytes'], url=m['acquisition_member']['url'], origin_plan_sha256=m['acquisition_origin_plan_sha256'],
    archived_munged_sha256=m['archived_munged_sha256'], archived_munged_bytes=m['archived_munged_bytes']) for m in new['members']]
receipt = dict(schema='independent_actual_extension7_metadata_closure_v1', all_checks_pass=True, fixed_current_metadata_sha256=fixed,
    source_member_count=100, commands=400, estimates=0, complete_scientific_members_unchanged_after_private_namespace_normalization=True,
    source_ledger=source_ledger, raw_declared_bytes=raw_bytes, archived_processed_declared_bytes=archived_bytes,
    scientific_source_version_origin_estimator_QC_guard_and_reservation_fields_unchanged=science_keys,
    all494_inherited_dependency_identities_unchanged_but_not_freshly_verified=True,
    all24_new_dependency_metadata_hashes_current=fresh, required_before_root_admission_full_current_dependency_and_archived_input_rehash=True,
    current_preserved_author_rejection_and_old_preflight_seals=history, unchanged_top_level_AST=unchanged,
    unchanged_nested_AST=['limits', 'regular_output_hashes', 'validation'], current_output_attempt_pending_seal_master_absent=absent,
    runtime_reference_source_body_reads=0, worker_estimator_network_mutex_operations=0, original1729_payloads_or402_controls_or190_reruns=0,
    execution_admission=False,
    qualification='Runner retains checkpoint operation and enforces full100 at final terminal. Frozen v7 execution preconditions additionally require root to verify genuine complete100 before initial admission; that root ordering condition is not supplied by this metadata review or by the runner initial gate.')
with (R / 'independent_extension_pipeline_v7_metadata_receipt_v1.json').open('x') as stream:
    json.dump(receipt, stream, indent=2, allow_nan=False)
    stream.write('\n')
print(json.dumps({'all_pass': True, 'members': 100, 'new_metadata_hashes': 24, 'inherited_deferred': 494}))

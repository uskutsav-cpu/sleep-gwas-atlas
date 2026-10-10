"""Metadata/AST draft closure only; never invokes prepare or execute."""
import ast
import hashlib
import json
import os
from pathlib import Path

P = Path(__file__).resolve().parents[1]
S = P / 'scripts'
R = P / 'reviews'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def tree(name):
    return ast.parse((S / name).read_text())

def functions(name):
    return {node.name: ast.dump(node, include_attributes=False) for node in tree(name).body if isinstance(node, ast.FunctionDef)}

def assignments(parsed, name):
    return [ast.dump(node.value, include_attributes=False) for node in ast.walk(parsed) if isinstance(node, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == name for t in node.targets)]

def h2update(parsed):
    return [ast.dump(node, include_attributes=False) for node in ast.walk(parsed) if isinstance(node, ast.Expr)
        and isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Attribute)
        and isinstance(node.value.func.value, ast.Name) and node.value.func.value.id == 'job' and node.value.func.attr == 'update']

fixed = {
    str(S / '58_run_finngen_feasibility_stage_v3_9.py'): 'b782d30e951fdde4c69c0c461e4f5e422332103446e6eb77a5d08ad4aa41710c',
    str(S / '59_prepare_finngen_feasibility_stage_v3_9.py'): '32b921dff9228f061d777a76fffb6190f9bdd19f6e899f00d21bb1e3ab9f8c39',
    str(S / '60_verify_finngen_h2_diagnostic_v5.py'): 'e9e85f9a79d6cc07da4da9b41b3a0d6bd9b95a969b4165d31a0217e0950719e2',
    str(P / 'manifests/finngen_preprocessing_plan_v4_3_8.json'): '78f0faea56024be8eca65421bf655c0f7c224afe6f029ed3ca0cea22491781b8',
    str(P / 'manifests/finngen_preprocessing_admission_v4_3_8.json'): '4c45bcdc31d23f9370c20dea7d6610894dc5fec57f229898a477befea039ab99',
}
assert all(sha(path) == digest for path, digest in fixed.items())
f6 = functions('58_run_finngen_feasibility_stage_v3_6.py')
f8 = functions('58_run_finngen_feasibility_stage_v3_8.py')
f9 = functions('58_run_finngen_feasibility_stage_v3_9.py')
old, new = tree('59_prepare_finngen_feasibility_stage_v3_6.py'), tree('59_prepare_finngen_feasibility_stage_v3_9.py')
checks = dict(
    exact_qualified3_8_worker_reader_AST=f8['read_fixed_worker'] == f9['read_fixed_worker'],
    exact_qualified3_8_execute_AST=f8['execute'] == f9['execute'],
    exact_prior_producer_gate3_6_AST=f6['prior_preprocessing_gate'] == f9['prior_preprocessing_gate'],
    same_h2_science_job_update_AST=h2update(old) == h2update(new),
    same_assumed_source_identity_and_effectiveN_AST=assignments(old, 'identity') == assignments(new, 'identity'),
    same_original_expected_source_AST=assignments(old, 'expected') == assignments(new, 'expected'),
    same_baseline_plan_and_SHA_AST=assignments(old, 'BASE') == assignments(new, 'BASE') and assignments(old, 'BASE_SHA') == assignments(new, 'BASE_SHA'),
)
text = (S / '59_prepare_finngen_feasibility_stage_v3_9.py').read_text()
checks['fixed_actual3_8_plan_SHA_before_prior_parse'] = "if sha(prior_plan)!='78f0faea56024be8eca65421bf655c0f7c224afe6f029ed3ca0cea22491781b8'" in text
checks['fixed_prior_stage_result_snapshot_before_parse_and_after'] = text.index('frozen_prior=') < text.index('prior = json.loads(prior_stage.read_text())') and text.count('check_hashes(frozen_prior)') == 2
checks['only_observed_h2_CLI'] = "choices=['observed_h2']" in text
checks['active_prior_paths_point3_8'] = all(x in text for x in [
    "prior_plan = P/'manifests/finngen_preprocessing_plan_v4_3_8.json'",
    "prior_stage = P/'logs/finngen_preprocessing_stage_receipt_v4_3_8.json'",
    "derivative = SOURCE_NAMESPACE/'pipeline_replay_v3_8/derived/insomnia.sumstats.gz'",
    "preprocessing_receipt = SOURCE_NAMESPACE/'pipeline_replay_v3_8/derived/insomnia.preprocessing.json'",
    "P/'scripts/59_prepare_finngen_feasibility_stage_v3_8.py'",
])
assert all(checks.values()), checks
seals = [
    ('independent_finngen_h2_candidate_findings_seal_v4_3_3.json', 'a7e236c85b86d58c40b8612f37a4fd1d32790d7729854ccbbd968873a05a4c8f', 474),
    ('independent_finngen_h2_corrected_candidate_findings_seal_v4_3_4.json', '717aaff9ee2cb7e646dcf4ab492e222388ff2149958efc16d4b6689327f32bc3', 405),
    ('independent_finngen_h2_report_commit_findings_seal_v4_3_5.json', '8542d2fd6d2a3228941317a4c7920e2d1adc212fe23ca2bbce3e1dcad3126c35', 296),
    ('independent_finngen_h2_final_consumer_review_seal_v4_3_6.json', 'facb1f7253b036dddc07d69ed26c484871cafb6cbdf2a14188f3bc4105785b3f', 62),
    ('independent_finngen_operational_prelaunch_seal_v4_3_8.json', '903b4d89868e19b6876311637ff7a8d37fc3383a57cded16abdf47cf96f8322f', 943),
]
history = []
for name, digest, count in seals:
    path = R / name
    assert sha(path) == digest
    record = json.loads(path.read_text())
    assert len(record['artifact_sha256']) == count
    assert all(sha(p) == h for p, h in record['artifact_sha256'].items())
    history.append(dict(seal=str(path), seal_sha256=digest, artifact_count=count, current_regular_artifact_hashes_match=True))
prospective = [P / 'manifests/finngen_observed_h2_plan_v4_3_9.json', P / 'manifests/finngen_observed_h2_admission_v4_3_9.json',
    P / 'logs/finngen_observed_h2_stage_receipt_v4_3_9.json']
absent = {str(path): not (path.exists() or path.is_symlink()) for path in prospective}
assert all(absent.values())
actual_preprocessing = json.loads((P / 'manifests/finngen_preprocessing_plan_v4_3_8.json').read_text())
source_fields = {k: actual_preprocessing[k] for k in ['source', 'source_identity', 'source_generation', 'assumed_effective_N', 'input_sha256', 'guard', 'global_reservation_bytes']}
receipt = dict(schema='independent_finngen_h2_prospective_closure_v4_3_9', fixed_candidate_producer_consumer_sha256=fixed,
    AST_and_binding_checks=checks, all_checks_pass=True, current_preserved_seals=history,
    old_rejected_artifact_count=1175, accepted_final_consumer_artifact_count=62, accepted_preprocessor_artifact_count=943,
    future_actual_producer_metadata_source=source_fields, biological_source_body_read=False,
    actual_h2_plan_admission_stage_absent=absent, current3_8_preprocessing_root_admission_verified=True,
    qualified_science_and_worker_controls_inherited=True, actual_h2_execution_admission=False,
    prospective_preparer_not_executed=True, fits_workers_locks_network=0)
with (R / 'independent_finngen_h2_candidate_closure_receipt_v4_3_9.json').open('x') as stream:
    json.dump(receipt, stream, indent=2, allow_nan=False)
    stream.write('\n')
print(json.dumps({'all_checks_pass': True, 'history_artifact_entries': sum(row['artifact_count'] for row in history), 'actual_h2_plan_absent': True}))

"""Bounded AST/command donor comparison; no scientific body reads or workers."""
from pathlib import Path
import ast
import difflib
import json

P = Path(__file__).resolve().parents[1]
R = P / 'reviews'
a = (P / 'scripts/58_run_finngen_feasibility_stage_v3_7.py').read_text()
b = (P / 'scripts/58_run_finngen_feasibility_stage_v3_8.py').read_text()

def norm(text):
    return text.replace('v4_3_7', 'v4_3_8').replace('v3_7', 'v3_8')

def functions(text):
    return {n.name: ast.dump(n, include_attributes=False) for n in ast.parse(text).body if isinstance(n, ast.FunctionDef)}

af, bf = functions(norm(a)), functions(b)
unchanged = [k for k in af if af[k] == bf[k]]
changed = [k for k in af if af[k] != bf[k]]
added = sorted(set(bf) - set(af))
assert changed == ['execute'] and added == ['read_fixed_worker']
oldexecute = next(n for n in ast.parse(norm(a)).body if isinstance(n, ast.FunctionDef) and n.name == 'execute')
newexecute = next(n for n in ast.parse(b).body if isinstance(n, ast.FunctionDef) and n.name == 'execute')
oldlimits = next(n for n in ast.walk(oldexecute) if isinstance(n, ast.FunctionDef) and n.name == 'limits')
newlimits = next(n for n in ast.walk(newexecute) if isinstance(n, ast.FunctionDef) and n.name == 'limits')
assert ast.dump(oldlimits, include_attributes=False) == ast.dump(newlimits, include_attributes=False)
oldplan = json.loads((P / 'manifests/finngen_preprocessing_plan_v4_3_7.json').read_text())
newplan = json.loads((P / 'manifests/finngen_preprocessing_plan_v4_3_8.json').read_text())
worker = str(P / 'scripts/56_preprocess_finngen_insomnia_feasibility_v3.py')
receipt = dict(
    schema='independent_finngen_exact_additive_code_correction_v4_3_8',
    namespace_normalized_top_level_AST_unchanged=unchanged,
    changed_functions=changed, added_functions=added, resource_meter_AST_unchanged=True,
    scientific_worker_SHA_unchanged=oldplan['dependencies_sha256'][worker] == newplan['dependencies_sha256'][worker],
    all170_donor_dependencies_retained=all(newplan['dependencies_sha256'].get(k) == v for k, v in oldplan['dependencies_sha256'].items()),
    all_science_argv_equal_except_plan_epoch=norm(json.dumps(oldplan['jobs'][0]['command_template'])) == json.dumps(newplan['jobs'][0]['command_template']),
    executor_diff='\n'.join(difflib.unified_diff(a.splitlines(), b.splitlines(), fromfile='58_3_7', tofile='58_3_8')),
    checker_correction='Initial unsaved inline probe used a relative worker key and failed KeyError before any artifact write. This saved checker uses exact absolute plan keys; no candidate defect or mutation occurred.',
    real_worker_fits_body_reads_network_lock=False)
assert receipt['scientific_worker_SHA_unchanged'] and receipt['all170_donor_dependencies_retained'] and receipt['all_science_argv_equal_except_plan_epoch']
with (R / 'independent_finngen_code_correction_receipt_v4_3_8.json').open('x') as f:
    json.dump(receipt, f, indent=2)
    f.write('\n')
print(json.dumps({'unchanged_functions': unchanged, 'changed_functions': changed, 'added_functions': added}))

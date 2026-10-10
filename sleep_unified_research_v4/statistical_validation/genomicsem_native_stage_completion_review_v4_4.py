#!/usr/bin/env python3
"""Narrow independent first-observation controls; synthetic metadata only."""
import ast
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import sys
import tempfile

P = Path(__file__).resolve().parents[1]
S = P / 'scripts'
V = P / 'statistical_validation'
R = P / 'reviews'
SELECTOR = S / 'native_stage_completion_v4_3.py'
OLD = S / '33_compare_native_campaign_v4_3.py'
NEW = S / '33_compare_native_campaign_v4_4.py'
EXPECTED = {str(SELECTOR): 'ed3b4d7af808d004447c955bf1ebafe6f70951ce9a178cab000ba865e13d487c',
            str(OLD): 'abf7799709309f86af415badb244f810145de511d91f99f231bdce0a0356c47f',
            str(NEW): '391bfa7863bbbca3c971523b52e1f87cda2b71e1a8da83b094ab2ffd576493e4',
            str(R / 'genomicsem_native_stage_completion_review_seal_v4_3.json'): '24d0a2ac3b7935dc49a37f9c8a18c7611f5d26faffffbfe2dcf84a9abdc9945a'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run():
    for path, expected in EXPECTED.items():
        assert sha(path) == expected, path
    old_seal = json.loads((R / 'genomicsem_native_stage_completion_review_seal_v4_3.json').read_text())
    for p, meta in old_seal['artifacts'].items():
        assert sha(p) == meta['sha256'] and Path(p).stat().st_size == meta['bytes'], p
    a, b = ast.parse(OLD.read_text()), ast.parse(NEW.read_text())
    a.body = [n for n in a.body if not isinstance(n, ast.FunctionDef) or n.name != 'record']
    b.body = [n for n in b.body if not isinstance(n, ast.FunctionDef) or n.name != 'record']
    assert ast.dump(a, include_attributes=False) == ast.dump(b, include_attributes=False)
    sys.path.insert(0, str(S))
    fixture = load_module(V / 'genomicsem_native_stage_completion_review_v4_3.py', '_original_stage_fixture_v4_4')
    mod = load_module(SELECTOR, '_selector_first_binding_v4_4')
    col = load_module(NEW, '_collator_first_binding_v4_4')
    checks = [{'case': 'entire_program_AST_unchanged_except_record', 'pass': True},
              {'case': 'original_four_artifacts_and_seal_unchanged', 'pass': True}]
    with tempfile.TemporaryDirectory(prefix='native-first-binding-review-') as td:
        root = Path(td)
        f = root / 'small_metadata.txt'
        f.write_text('AAAA')
        col.INPUTS = {}
        assert col.record(f) == f
        original = dict(col.INPUTS[str(f)])
        original_object = col.INPUTS[str(f)]
        assert col.record(f) == f and col.INPUTS[str(f)] is original_object
        checks.append({'case': 'identical_repeat_retains_first_object_hash_and_size', 'pass': True})
        f.write_text('BBBB')
        try:
            col.record(f)
        except RuntimeError as exc:
            assert 'CONSUMED_SOURCE_CHANGED_SINCE_FIRST_OBSERVATION' in str(exc)
        else:
            raise AssertionError('SAME_SIZE_DRIFT_NOT_REJECTED')
        assert col.INPUTS[str(f)] == original
        checks.append({'case': 'same_size_changed_hash_rejected_and_first_identity_retained', 'pass': True})
        f.write_text('BBBBB')
        try:
            col.record(f)
        except RuntimeError as exc:
            assert 'CONSUMED_SOURCE_CHANGED_SINCE_FIRST_OBSERVATION' in str(exc)
        else:
            raise AssertionError('SIZE_DRIFT_NOT_REJECTED')
        assert col.INPUTS[str(f)] == original
        checks.append({'case': 'changed_size_rejected_and_first_identity_retained', 'pass': True})
        for target in ('selected_monitor', 'native_job_receipt'):
            f = fixture.fixture(mod, root / target)
            col.INPUTS = {}
            mod.stage_monitor('extension', record=col.record)
            before = {p: dict(v) for p, v in col.INPUTS.items()}
            if target == 'native_job_receipt':
                victim = Path(next(iter(f['mon']['native_receipt_sha256'])))
                r = json.loads(victim.read_text())
                r['audit_note'] = 'changed between gate calls'
                f['mon']['native_receipt_sha256'][str(victim)] = fixture.dump(victim, r)
            else:
                f['mon']['audit_note'] = 'changed between gate calls'
            fixture.dump(f['mon_path'], f['mon'])
            try:
                mod.stage_monitor('extension', record=col.record)
            except RuntimeError as exc:
                assert 'CONSUMED_SOURCE_CHANGED_SINCE_FIRST_OBSERVATION' in str(exc)
            else:
                raise AssertionError('ORIGINAL_WITNESS_NOT_REJECTED: ' + target)
            assert col.INPUTS == before
            assert not all(col.sha(p) == m['sha256'] for p, m in col.INPUTS.items())
            checks.append({'case': 'original_' + target + '_witness_now_rejected_before_output', 'pass': True})
        f = fixture.fixture(mod, root / 'stable_complete_metadata')
        col.INPUTS = {}
        mod.stage_monitor('extension', record=col.record)
        before = {p: dict(v) for p, v in col.INPUTS.items()}
        mod.stage_monitor('extension', record=col.record)
        assert col.INPUTS == before
        assert all(col.sha(p) == m['sha256'] for p, m in col.INPUTS.items())
        checks.append({'case': 'unchanged_two_gate_metadata_accepted_without_rebinding', 'pass': True})
    for path, expected in EXPECTED.items():
        assert sha(path) == expected, path
    for p, meta in old_seal['artifacts'].items():
        assert sha(p) == meta['sha256'], p
    return {'recorded_utc': datetime.now(timezone.utc).isoformat(),
            'status': 'NARROW_V4_4_FIRST_OBSERVATION_CORRECTION_PASS',
            'source_sha256_before_and_after': EXPECTED,
            'helper_sha256': sha(__file__), 'checks': checks, 'passing_control_count': len(checks),
            'prior_control_count_reused_by_exact_selector_and_seal_identity': 36,
            'prior_blocker_witness_count_rejected': 2,
            'scientific_arithmetic_tolerances_families_and_selector_unchanged': True,
            'original_rejection_preserved': True, 'collator_main_called': False,
            'real_workers_launched': 0, 'real_stage_monitor_consumed': False,
            'GWAS_body_or_genotype_reads': 0, 'actual_extension_completion_claim': False,
            'peak_review_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}


if __name__ == '__main__':
    out = V / 'genomicsem_native_stage_completion_review_controls_v4_4.json'
    if out.exists():
        raise RuntimeError('DISTINCT_REVIEW_OUTPUT_REQUIRED')
    result = run()
    with out.open('x') as f:
        f.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('status', 'passing_control_count', 'real_workers_launched')}))

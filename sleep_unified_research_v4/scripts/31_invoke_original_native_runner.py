#!/usr/bin/env python3
"""Invoke unchanged v1 estimator executor with an explicitly relocated PACKAGE."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
P = ROOT / 'sleep_unified_research_v4'
PLAN = P / 'manifests/ssd_native_execution_plan_v4_3.json'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    plan = json.loads(PLAN.read_text())
    original = ROOT / 'sleep_unified_research_v1/scripts/04_native_reproduction_runner.py'
    assert sha(original) == plan['original_runner_sha256'] == 'fa5ca8528caa1037cca55d5234e7e7721582abce479c95949bb81255e3428fa9'
    assert sha(Path(__file__)) == plan['invoker_sha256']
    for path, expected in plan['support_file_sha256'].items():
        assert sha(Path(path)) == expected, path
    spec = importlib.util.spec_from_file_location('unchanged_sleep_native_v1', original)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # Sole operational relocation. ROOT, inputs, statistics and preflight are unchanged.
    module.PACKAGE = Path(plan['ssd_support_package'])
    actual = module.jobs()
    assert actual == plan['jobs'], 'FROZEN_JOB_IDENTITY_CHANGED'
    if len(sys.argv) != 2 or sys.argv[1] not in ('core', 'extension', 'validation'):
        raise SystemExit('REQUIRE_ONE_STAGE')
    stage = sys.argv[1]
    sys.argv = [str(original), '--stage', stage, '--execute']
    module.main()

if __name__ == '__main__':
    main()

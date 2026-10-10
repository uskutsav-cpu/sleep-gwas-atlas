#!/usr/bin/env python3
"""Seal distinct preparation artifacts and verify prior specialist seal unchanged."""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    root = Path(__file__).resolve().parents[2]
    package = root / 'sleep_unified_research_v4'
    old = package / 'reviews/genomicsem_specialist_review_receipt_v4.json'
    previous = json.loads(old.read_text())
    # Preserve the previously sealed scientific review; inspect its schema
    # explicitly rather than silently accepting a failed historical recheck.
    artifacts = previous['records']
    for item in artifacts:
        path = Path(item['path'])
        if not path.is_absolute():
            path = root / path
        if sha(path) != item['sha256']:
            raise RuntimeError('PREVIOUS_SPECIALIST_SEAL_CHANGED: ' + str(path))
    paths = list((package / 'scripts').glob('canonical_*.py')) + [package / 'scripts/42_prepare_native_canonical_calibration.py',
            package / 'manifests/native_canonical_calibration_plan_v4.json', package / 'manifests/native_canonical_calibration_admission_v4.json',
            package / 'statistical_validation/canonical_native_calibration_fixture_receipt_v4.json',
            package / 'statistical_validation/canonical_native_calibration_fixture_receipt_v4_2.json',
            package / 'statistical_validation/canonical_native_calibration_admission_guard_receipt_v4.json',
            package / 'reviews/genomicsem_native_calibration_executor_review_v4.md', Path(__file__)]
    current = [{'path': str(p), 'bytes': p.stat().st_size, 'sha256': sha(p)} for p in sorted(paths)]
    record = {'schema': 'native_canonical_calibration_preparation_seal_v1', 'sealed_utc': datetime.now(timezone.utc).isoformat(),
              'artifacts': current, 'previous_specialist_review_receipt_sha256': sha(old), 'previous_specialist_seal_unchanged': True,
              'upstream_builder_receipt_sha256': sha(package / 'manifests/canonical_200_interval_preparation_v4.json'),
              'status': 'PREPARATION_AND_MOCK_CHECKS_ONLY_EXECUTION_BLOCKED', 'native_jobs_launched': 0,
              'independent_executor_binding_review_pass': False, 'realistic_LD_sampling_calibration_pass': False,
              'allow_41_covariance_outcomes': False, 'calibrated_biological_p_values_computed': False}
    path = package / 'reviews/genomicsem_native_calibration_preparation_receipt_v4.json'
    with path.open('x') as f:
        json.dump(record, f, indent=2, allow_nan=False)
        f.write('\n')
    print(json.dumps({'artifacts_sealed': len(current), 'bytes': sum(r['bytes'] for r in current), 'receipt_sha256': sha(path),
                      'previous_specialist_seal_unchanged': True, 'native_jobs_launched': 0}))


if __name__ == '__main__':
    main()

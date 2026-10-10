#!/usr/bin/env python3
"""Seal third preparation without editing any preceding sealed artifact."""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024**2), b''): h.update(block)
    return h.hexdigest()


def main():
    root = Path(__file__).resolve().parents[2]
    package = root / 'sleep_unified_research_v4'
    previous = [package / 'reviews/genomicsem_specialist_review_receipt_v4.json',
                package / 'reviews/genomicsem_native_calibration_preparation_receipt_v4.json',
                package / 'reviews/genomicsem_native_calibration_preparation_receipt_v4_2.json']
    for receipt in previous:
        record = json.loads(receipt.read_text())
        for item in record.get('records', record.get('artifacts', [])):
            path = Path(item['path'])
            if not path.is_absolute(): path = root / path
            if sha(path) != item['sha256']: raise RuntimeError('PREVIOUS_SEALED_ARTIFACT_CHANGED: ' + str(path))
    paths = list((package / 'scripts').glob('canonical_*_v4_3.py')) + [package / 'scripts/42_prepare_native_canonical_calibration_v4_3.py',
             package / 'manifests/native_canonical_calibration_plan_v4_3.json', package / 'manifests/native_canonical_calibration_admission_v4_3.json',
             package / 'statistical_validation/canonical_native_calibration_fixture_receipt_v4_4.json',
             package / 'statistical_validation/canonical_operational_correction_fixture_receipt_v4_3.json',
             package / 'statistical_validation/canonical_cleanup_and_final_guard_fixture_receipt_v4_3.json',
             package / 'statistical_validation/canonical_signal_diagnostic_post_hash_fixture_receipt_v4_3.json',
             package / 'statistical_validation/canonical_native_calibration_admission_guard_receipt_v4_3.json',
             package / 'reviews/genomicsem_native_calibration_executor_review_v4_3.md', Path(__file__)]
    artifacts = [{'path': str(p), 'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted(paths)]
    record = {'schema': 'native_canonical_calibration_third_preparation_seal_v1', 'sealed_utc': datetime.now(timezone.utc).isoformat(),
              'artifacts': artifacts, 'previous_receipt_sha256': {str(p): sha(p) for p in previous}, 'all_previous_seals_unchanged': True,
              'status': 'THIRD_OPERATIONAL_PREPARATION_EXECUTION_BLOCKED', 'native_jobs_launched': 0,
              'independent_third_executor_binding_review_pass': False, 'realistic_LD_sampling_calibration_pass': False,
              'allow_41_covariance_outcomes': False, 'calibrated_biological_p_values_computed': False}
    path = package / 'reviews/genomicsem_native_calibration_preparation_receipt_v4_3.json'
    with path.open('x') as f:
        json.dump(record, f, indent=2, allow_nan=False)
        f.write('\n')
    print(json.dumps({'artifacts_sealed': len(artifacts), 'bytes': sum(r['bytes'] for r in artifacts), 'receipt_sha256': sha(path),
                      'all_previous_seals_unchanged': True, 'native_jobs_launched': 0}))


if __name__ == '__main__': main()

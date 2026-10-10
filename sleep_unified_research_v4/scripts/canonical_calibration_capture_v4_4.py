#!/usr/bin/env python3
"""One future admitted calibration fit, with exact weighted-row capture."""
import argparse
import json
from pathlib import Path
import runpy
import sys

from canonical_calibration_common_v4_4 import admit, check_hashes, clean, sha, utc, validate_inherited_lock, write_new


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--admission', type=Path, required=True)
    p.add_argument('--plan-sha', required=True)
    p.add_argument('--job-id', required=True)
    p.add_argument('--shared-lock-fd', type=int, required=True)
    args = p.parse_args()
    validate_inherited_lock(args.shared_lock_fd)
    if sha(args.plan) != args.plan_sha:
        raise RuntimeError('PREPARED_PLAN_CHANGED')
    plan, admission = admit(args.plan, args.admission)
    admission_before = sha(args.admission)
    job = next(j for j in plan['jobs'] if j['job_id'] == args.job_id)
    if job['control_id'] not in plan['allowed_control_ids'] or len(job['inputs']) != 2:
        raise RuntimeError('ONLY_TWO_PRESPECIFIED_CONTROL_PAIRS_ADMITTED')
    before = check_hashes(job['input_sha256'])
    dependencies_before = check_hashes(plan['dependencies_sha256'])
    code = Path(plan['ldsc_dir'])
    sys.path.insert(0, str(code))
    import numpy as np
    import pandas as pd
    import scipy
    import ldscore.sumstats as ss
    import ldscore.regressions as reg
    import ldscore.jackknife as jk
    from canonical_native_adapter import NativeCalibrationAdapter
    if not sys.version.startswith('3.9.23') or (np.__version__, pd.__version__, scipy.__version__) != ('1.21.5', '1.3.3', '1.7.3'):
        raise RuntimeError('PINNED_NUMERICAL_ENVIRONMENT_CHANGED')
    blocks = json.loads(Path(plan['block_receipt_path']).read_text())
    capture = NativeCalibrationAdapter(blocks, job['partition'], Path(job['output_dir']) / 'weighted_capture', np, admitted_scope=admission['scope'])
    old_estimate = ss.estimate_rg
    original_pair = ss._rg
    recorded = []

    def guarded_pair(rows, effective_args, log, M, refs, weight_name, index):
        expected = None if job['method'] == 'one_step' else 30
        if effective_args.two_step != expected or index != 0 or effective_args.no_intercept:
            raise RuntimeError('EFFECTIVE_METHOD_DIFFERS_FROM_FROZEN_ARM')
        if effective_args.intercept_h2 != [None, None] or effective_args.intercept_gencov != [None, None]:
            raise RuntimeError('FROZEN_FREE_INTERCEPTS_CHANGED')
        return original_pair(rows, effective_args, log, M, refs, weight_name, index)

    def estimate(arguments, log):
        result = old_estimate(arguments, log)
        if len(result) != 1 or result[0] is None:
            raise RuntimeError('SINGLE_CONTROL_ESTIMATOR_RETURN_REQUIRED')
        obj = result[0]
        recorded.append({'arguments': clean(vars(arguments), np),
                         'point_estimates': capture.attributes(obj, ('rg_ratio', 'rg_jknife', 'rg_se', '_negative_hsq')),
                         'effective_method_asserted_before_regression': True, 'free_intercepts_asserted_before_regression': True,
                         **capture.record()})
        return result

    ss.estimate_rg = estimate
    ss._rg = guarded_pair
    try:
        with capture.install(ss, reg, jk):
            sys.argv = [str(code / 'ldsc.py')] + job['ldsc_args']
            runpy.run_path(str(code / 'ldsc.py'), run_name='__main__')
    finally:
        ss.estimate_rg = old_estimate
        ss._rg = original_pair
    if len(recorded) != 1:
        raise RuntimeError('EXPECTED_EXACTLY_ONE_COMPLETED_CONTROL')
    after = check_hashes(job['input_sha256'])
    dependencies_after = check_hashes(plan['dependencies_sha256'])
    if before != after or dependencies_before != dependencies_after or sha(args.plan) != args.plan_sha or sha(args.admission) != admission_before:
        raise RuntimeError('SOURCE_OR_CODE_CHANGED_DURING_CONTROL')
    # This receipt is created only after stock main completes. Partial arrays and
    # the stock log remain inspectable after a failed constructor or native main.
    write_new(Path(job['output_dir']) / 'native_capture_receipt.json',
              {'schema': 'canonical_native_method_calibration_capture_v1', 'completed_utc': utc(), 'job': job,
               'plan_sha256': args.plan_sha, 'admission_sha256': admission_before, 'python': sys.version,
               'libraries': {'numpy': np.__version__, 'pandas': pd.__version__, 'scipy': scipy.__version__},
               'input_sha256_before': before, 'input_sha256_after': after,
               'dependency_sha256_before': dependencies_before, 'dependency_sha256_after': dependencies_after,
               'capture': recorded[0], 'scope': 'METHOD_CALIBRATION_ONLY', 'calibrated_biological_p_values_computed': False})


if __name__ == '__main__':
    main()

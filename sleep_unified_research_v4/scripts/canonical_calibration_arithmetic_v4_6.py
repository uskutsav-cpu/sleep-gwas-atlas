#!/usr/bin/env python3
"""Independent captured-row solve and duplicate/ratio arithmetic, no LDSC import.

Run only after all prespecified native controls complete. Uses raw retained rows
with np.linalg.lstsq, never LDSC block cross-product/deletion helper functions.
Weights and Nbar are the captured full-fit quantities, as in stock fast LDSC.
This is an implementation check, not realistic-LD sampling calibration.
"""
import argparse
import json
from pathlib import Path

from canonical_calibration_common_v4_6 import admit, sha, utc, validate_inherited_lock, write_new


def validate(plan, np):
    checks = []
    captures = {}
    point_atol, point_rtol = 1e-10, 1e-8
    cov_atol, cov_rtol = 1e-14, 1e-6

    def compare(name, actual, expected, covariance=False):
        a, b = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
        if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
            raise RuntimeError('ARITHMETIC_SHAPE_OR_FINITE_FAILURE: ' + name)
        tol = (cov_atol + cov_rtol * np.abs(b)) if covariance else (point_atol + point_rtol * np.abs(b))
        maximum = float(np.max(np.abs(a - b)))
        passed = bool(np.all(np.abs(a - b) <= tol))
        checks.append({'name': name, 'max_absolute_difference': maximum, 'pass': passed})
        if not passed:
            raise RuntimeError('INDEPENDENT_ARITHMETIC_DISAGREES: ' + name)

    def exact(name, left, right):
        passed = left == right
        checks.append({'name': name, 'exact_equal': bool(passed), 'pass': bool(passed)})
        if not passed:
            raise RuntimeError('CROSS_ARM_EXACT_IDENTITY_DISAGREES: ' + name)

    def load(item):
        if sha(item['path']) != item['sha256']:
            raise RuntimeError('CAPTURE_ARRAY_CHANGED')
        a = np.load(item['path'], allow_pickle=False)
        if list(a.shape) != item['shape'] or str(a.dtype) != item['dtype'] or not np.isfinite(a).all():
            raise RuntimeError('CAPTURE_ARRAY_INVALID')
        return a

    def covariance(point, deletes):
        pseudo = 200 * point - 199 * deletes
        centered = pseudo - pseudo.mean(axis=0)
        return centered.T @ centered / (200 * 199)

    for job in plan['jobs']:
        path = Path(job['output_dir']) / 'native_capture_receipt.json'
        r = json.loads(path.read_text())
        if r['job'] != job or r['plan_sha256'] != sha(plan['self_path']):
            raise RuntimeError('CONTROL_RECEIPT_IDENTITY_MISMATCH')
        capture = r['capture']
        captures[job['job_id']] = capture
        # Independently solve only replicate A; replicate B remains an exact
        # same-source implementation identity control against all captured data.
        if job['duplicate'] != 'A':
            continue
        totals, total_deletes = [], []
        for base in capture['regressions']:
            stages = []
            for stage in base['stages']:
                x, y = load(stage['final_weighted_design']), load(stage['final_weighted_response'])
                cuts = stage['actual_separators']
                if len(cuts) != 201 or cuts[0] != 0 or cuts[-1] != len(x) or any(a >= b for a, b in zip(cuts, cuts[1:])):
                    raise RuntimeError('NONEMPTY_BLOCK_REQUIREMENT_FAILED')
                beta, _, rank, _ = np.linalg.lstsq(x, y, rcond=None)
                if rank != x.shape[1]:
                    raise RuntimeError('FULL_WEIGHTED_DESIGN_RANK_DEFICIENT')
                deleted = []
                for left, right in zip(cuts, cuts[1:]):
                    keep_x = np.concatenate((x[:left], x[right:]), axis=0)
                    keep_y = np.concatenate((y[:left], y[right:]), axis=0)
                    b, _, rank, _ = np.linalg.lstsq(keep_x, keep_y, rcond=None)
                    if rank != x.shape[1]:
                        raise RuntimeError('DELETE_WEIGHTED_DESIGN_RANK_DEFICIENT')
                    deleted.append(b[:, 0])
                point, deletes = beta.T, np.asarray(deleted)
                prefix = job['job_id'] + '/' + base['name'] + '/' + stage['name']
                compare(prefix + '/row_point', point, stage['jackknife']['est'])
                compare(prefix + '/row_deletes', deletes, stage['jackknife']['delete_values'])
                compare(prefix + '/row_covariance', covariance(point, deletes), stage['jackknife']['jknife_cov'], covariance=True)
                stages.append((point, deletes))
            if base['two_step']:
                first, second = stages
                N = load(capture['pair']['N1']) if base['name'] == 'hsq1' else load(capture['pair']['N2'])
                if base['name'] == 'gencov':
                    N = np.sqrt(load(capture['pair']['N1']) * load(capture['pair']['N2']))
                xraw = load(capture['pair']['raw_reference_design'])
                initial = load(base['initial_weights'])
                scaled = N.reshape(-1, 1) * xraw / base['Nbar']
                c = np.sum(initial * scaled) / np.sum(initial * scaled**2)
                compare(job['job_id'] + '/' + base['name'] + '/independent_c', c, base['two_step_combination']['c'])
                intercept = first[0][0, -1]
                point = np.column_stack((second[0], np.array([[intercept]])))
                slopes = second[1] - c * (first[1][:, -1] - intercept).reshape(-1, 1)
                deletes = np.column_stack((slopes, first[1][:, -1]))
            else:
                point, deletes = stages[0]
            prefix = job['job_id'] + '/' + base['name']
            compare(prefix + '/combined_point', point, base['combined_jackknife']['est'])
            compare(prefix + '/combined_deletes', deletes, base['combined_jackknife']['delete_values'])
            compare(prefix + '/combined_covariance', covariance(point, deletes), base['combined_jackknife']['jknife_cov'], covariance=True)
            M = np.asarray(base['M'])
            Ncheck = load(capture['pair']['N1']) if base['name'] == 'hsq1' else load(capture['pair']['N2'])
            if base['name'] == 'gencov':
                Ncheck = np.sqrt(load(capture['pair']['N1']) * load(capture['pair']['N2']))
            compare(prefix + '/independent_Nbar', np.mean(Ncheck), base['Nbar'])
            compare(prefix + '/independent_M_tot', np.sum(M), base['M_tot'])
            total = point[:, :1] @ M.T / base['Nbar']
            deletion = deletes[:, :1] @ M.T / base['Nbar']
            compare(prefix + '/total', total.squeeze(), base['fitted_attributes']['tot'])
            compare(prefix + '/total_deletes', deletion, load(base['tot_delete_values']))
            totals.append(total.item())
            total_deletes.append(deletion)
        a, b, c = totals
        if a <= 0 or b <= 0 or np.any(total_deletes[0] <= 0) or np.any(total_deletes[1] <= 0):
            raise RuntimeError('NONPOSITIVE_HERITABILITY_OR_DELETE_RATIO_DENOMINATOR')
        ratio = np.array([[c / np.sqrt(a * b)]])
        ratio_deletes = total_deletes[2] / np.sqrt(total_deletes[0] * total_deletes[1])
        compare(job['job_id'] + '/independent_ratio', ratio, capture['ratio']['point_ratio'])
        compare(job['job_id'] + '/independent_ratio_cov', covariance(ratio, ratio_deletes), capture['ratio']['native_jackknife']['jknife_cov'], covariance=True)
        pseudo = 200 * ratio - 199 * ratio_deletes
        compare(job['job_id'] + '/independent_ratio_jknife', pseudo.mean(axis=0).reshape(1, 1), capture['ratio']['native_jackknife']['jknife_est'])

    for control in plan['allowed_control_ids']:
        for method in ('one_step', 'two_step30'):
            for partition in ('stock_default200', 'canonical200'):
                jobs = [j for j in plan['jobs'] if (j['control_id'], j['method'], j['partition']) == (control, method, partition)]
                left, right = [captures[j['job_id']] for j in jobs]
                prefix = control + '/' + method + '/' + partition
                if left['pair']['ordered_final_snp_sha256'] != right['pair']['ordered_final_snp_sha256']:
                    raise RuntimeError('DUPLICATE_FINAL_SNP_ORDER_CHANGED')
                for key in ('rg_ratio', 'rg_jknife', 'rg_se'):
                    compare(prefix + '/duplicate_' + key, left['point_estimates'][key], right['point_estimates'][key])
                for lb, rb in zip(left['regressions'], right['regressions']):
                    compare(prefix + '/' + lb['name'] + '/duplicate_total_deletes', load(lb['tot_delete_values']), load(rb['tot_delete_values']))
                    if [s['actual_separators'] for s in lb['stages']] != [s['actual_separators'] for s in rb['stages']]:
                        raise RuntimeError('DUPLICATE_BLOCK_IDENTITY_CHANGED')
                    for ls, rs in zip(lb['stages'], rb['stages']):
                        if ls['final_weighted_design']['sha256'] != rs['final_weighted_design']['sha256'] or ls['final_weighted_response']['sha256'] != rs['final_weighted_response']['sha256']:
                            raise RuntimeError('DUPLICATE_WEIGHTED_ROWS_NOT_EXACTLY_IDENTICAL')
                def ratio_deletes(cap):
                    return load(cap['ratio']['numerator_delete']) / load(cap['ratio']['denominator_delete'])
                dl, dr = ratio_deletes(left), ratio_deletes(right)
                point = np.array([[left['point_estimates']['rg_ratio'], right['point_estimates']['rg_ratio']]])
                paired_cov = covariance(point, np.column_stack((dl, dr)))
                compare(prefix + '/duplicate_cov_equals_variance', paired_cov, np.full((2, 2), paired_cov[0, 0]), covariance=True)
                delta = float(point[0, 0] - point[0, 1])
                direct_variance = covariance(np.array([[delta]]), dl - dr).item()
                subtraction_variance = float(paired_cov[0, 0] + paired_cov[1, 1] - 2 * paired_cov[0, 1])
                if abs(delta) > 1e-12 or abs(direct_variance) > 1e-24 or abs(subtraction_variance) > 1e-24:
                    raise RuntimeError('KNOWN_ZERO_CONTRAST_IDENTITY_FAILED')
                checks.append({'name': prefix + '/known_zero_contrast', 'delta': delta, 'direct_delete_delta_variance': direct_variance,
                               'covariance_subtraction_variance': subtraction_variance, 'se': float(np.sqrt(direct_variance)), 'p_value': None, 'reason': '0/0 identity; no hypothesis P value', 'pass': True})
            default = next(captures[j['job_id']] for j in plan['jobs'] if (j['control_id'], j['method'], j['partition'], j['duplicate']) == (control, method, 'stock_default200', 'A'))
            canonical = next(captures[j['job_id']] for j in plan['jobs'] if (j['control_id'], j['method'], j['partition'], j['duplicate']) == (control, method, 'canonical200', 'A'))
            prefix = control + '/' + method + '/cross_arm'
            exact(prefix + '/final_ordered_SNP_sha256', canonical['pair']['ordered_final_snp_sha256'], default['pair']['ordered_final_snp_sha256'])
            exact(prefix + '/final_SNP_count', canonical['pair']['ordered_final_snp_count'], default['pair']['ordered_final_snp_count'])
            exact(prefix + '/effective_two_step', canonical['pair']['args_two_step'], default['pair']['args_two_step'])
            for key in ('N1', 'N2', 'raw_reference_design', 'canonical_block_labels'):
                exact(prefix + '/' + key + '_sha256', canonical['pair'][key]['sha256'], default['pair'][key]['sha256'])
            compare(control + '/' + method + '/default_vs_canonical_ratio_point', canonical['point_estimates']['rg_ratio'], default['point_estimates']['rg_ratio'])
            for db, cb in zip(default['regressions'], canonical['regressions']):
                bp = prefix + '/' + db['name']
                for key in ('name', 'Nbar', 'M', 'M_tot', 'two_step'):
                    exact(bp + '/' + key, cb[key], db[key])
                if cb['two_step']:
                    exact(bp + '/first_stage_ordered_row_indices_sha256', cb['first_stage_row_indices']['sha256'], db['first_stage_row_indices']['sha256'])
                checks.append({'name': bp + '/initial_weights/hash_observation', 'canonical_sha256': cb['initial_weights']['sha256'],
                               'default_sha256': db['initial_weights']['sha256'], 'binary_equal': cb['initial_weights']['sha256'] == db['initial_weights']['sha256'],
                               'pass': True, 'qualification': 'Native Gencov weights depend on previously fitted h2 points/intercepts and may inherit binary64 block-sum roundoff.'})
                compare(bp + '/initial_weights/raw_array_numerical_identity', load(cb['initial_weights']), load(db['initial_weights']))
                exact(bp + '/stage_count', len(cb['stages']), len(db['stages']))
                for cs, ds in zip(cb['stages'], db['stages']):
                    sp = bp + '/' + cs['name']
                    exact(sp + '/stage_name', cs['name'], ds['name'])
                    exact(sp + '/canonical_mask_block_counts', cs['canonical_block_counts'], ds['canonical_block_counts'])
                    for key in ('final_weighted_design', 'final_weighted_response'):
                        binary_equal = cs[key]['sha256'] == ds[key]['sha256']
                        checks.append({'name': sp + '/' + key + '/hash_observation', 'canonical_sha256': cs[key]['sha256'],
                                       'default_sha256': ds[key]['sha256'], 'binary_equal': binary_equal, 'pass': True,
                                       'qualification': 'Native sequential Gencov and second-stage response/weights may inherit binary64 block-sum roundoff; unequal hashes require numerical check under the frozen point/delete tolerance.'})
                        # Load/validate both arrays even when identities match.
                        # This keeps captured rows bound to the comparison.
                        compare(sp + '/' + key + '/raw_array_numerical_identity', load(cs[key]), load(ds[key]))
                compare(control + '/' + method + '/' + db['name'] + '/default_vs_canonical_total_point', cb['fitted_attributes']['tot'], db['fitted_attributes']['tot'])
                compare(control + '/' + method + '/' + db['name'] + '/default_vs_canonical_intercept_point', cb['fitted_attributes']['intercept'], db['fitted_attributes']['intercept'])
    return {'schema': 'canonical_independent_control_arithmetic_v1', 'completed_utc': utc(), 'checks': checks, 'all_checks_pass': all(c['pass'] for c in checks),
            'scope': 'METHOD_CALIBRATION_IMPLEMENTATION_CHECK_ONLY', 'realistic_LD_sampling_calibration_pass': False,
            'allow_41_covariance_outcomes': False, 'calibrated_biological_p_values_computed': False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--admission', type=Path, required=True)
    parser.add_argument('--shared-lock-fd', type=int, required=True)
    args = parser.parse_args()
    validate_inherited_lock(args.shared_lock_fd)
    plan, _ = admit(args.plan, args.admission)
    import numpy as np
    write_new(Path(plan['ssd_output_root']) / 'independent_control_arithmetic.json', validate(plan, np))


if __name__ == '__main__':
    main()

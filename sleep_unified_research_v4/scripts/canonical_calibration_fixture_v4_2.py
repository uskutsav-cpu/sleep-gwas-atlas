#!/usr/bin/env python3
"""Tiny synthetic adapter/arithmetic fixtures. No LDSC numerical module import.

Uses injected mock estimator interfaces and 600 invented SNP IDs. Fixture output
cannot admit execution or realistic-LD uncertainty calibration.
"""
import csv
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace

from canonical_calibration_common_v4_2 import PACKAGE, sha, utc, write_new


def main():
    import numpy as np
    import pandas as pd
    from canonical_native_adapter import NativeCalibrationAdapter, prototype
    from canonical_calibration_arithmetic_v4_2 import validate
    module = prototype()
    interfaces = module.verify_pinned_interfaces()
    events = []
    ids = ['fixture_rs_%d_%d' % (b, j) for b in range(200) for j in (1, 2, 3)]
    coordinates = {s: (1, 10 * (i // 3) + i % 3 + 1) for i, s in enumerate(ids)}
    rows = pd.DataFrame({'SNP': ids, 'N1': np.full(600, 100000.), 'N2': np.full(600, 90000.),
                         'Z1': [np.sqrt((40 if j == 1 else 2) + .01*b) for b in range(200) for j in (1, 2, 3)],
                         'Z2': [np.sqrt((40 if j == 2 else 2) + .015*b) for b in range(200) for j in (1, 2, 3)],
                         'L2': [2 + .01*b + .0001*j for b in range(200) for j in (1, 2, 3)], 'W_L2': np.ones(600)})
    reg, jk = SimpleNamespace(), SimpleNamespace()

    def jknife(point, deletes):
        pseudo = 200 * point - 199 * deletes
        cov = np.atleast_2d(np.cov(pseudo.T, ddof=1) / 200)
        return {'est': point, 'delete_values': deletes, 'jknife_est': pseudo.mean(axis=0).reshape(1, -1),
                'jknife_cov': cov, 'jknife_var': np.diag(cov).reshape(1, -1), 'jknife_se': np.sqrt(np.diag(cov)).reshape(1, -1)}

    class MockFast:
        def __init__(self, x, y, n_blocks=None, separators=None):
            cuts = np.floor(np.linspace(0, len(x), n_blocks + 1)).astype(int).tolist() if separators is None else list(separators)
            self.separators = cuts
            # Mock fast block cross-products deliberately differ from the
            # independent raw-row lstsq verification implementation.
            xtx, xty = x.T @ x, x.T @ y
            point = np.linalg.solve(xtx, xty).T
            deletes = np.asarray([np.linalg.solve(xtx - x[l:r].T @ x[l:r], xty - x[l:r].T @ y[l:r])[:, 0] for l, r in zip(cuts, cuts[1:])])
            self.__dict__.update(jknife(point, deletes))

    class MockRatio:
        def __init__(self, est, numer_delete_values, denom_delete_values):
            self.__dict__.update(jknife(est, numer_delete_values / denom_delete_values))

    class MockIRWLS:
        def __init__(self, x, y, update_func, n_blocks, w=None, slow=False, separators=None):
            result = jk.LstsqJackknifeFast(x * np.sqrt(w), y * np.sqrt(w), n_blocks=n_blocks, separators=separators)
            self.__dict__.update(result.__dict__)

    class MockBase:
        def __init__(self, y, x, w, N, M, n_blocks, intercept=None, slow=False, step1_ii=None, old_weights=False):
            initial = self._update_weights(x, w, N, float(M.sum()), None, intercept)
            Nbar = np.mean(N)
            scaled = N * x / Nbar
            design = np.column_stack((scaled, np.ones(len(scaled))))
            if step1_ii is None:
                fit = reg.IRWLS(design, y, None, n_blocks, w=initial)
                self.twostep_filtered = None
            else:
                mask = step1_ii.reshape(-1)
                first = reg.IRWLS(design[mask], y[mask], None, n_blocks, w=initial[mask])
                mapping = np.flatnonzero(mask)
                incoming = [0] + [int(mapping[i]) for i in first.separators[1:-1]] + [len(mask)]
                second = reg.IRWLS(scaled, y - first.est[0, -1], None, n_blocks, w=initial, separators=incoming)
                c = np.sum(initial * scaled) / np.sum(initial * scaled**2)
                fit = self._combine_twostep_jknives(first, second, float(M.sum()), c, Nbar)
                self.twostep_filtered = len(mask) - int(mask.sum())
            self.jknife = fit
            self.tot = float(fit.est[0, 0] * M[0, 0] / Nbar)
            self.tot_cov = float(fit.jknife_cov[0, 0] * (M[0, 0] / Nbar)**2)
            self.tot_se = float(np.sqrt(self.tot_cov))
            self.coef, self.coef_se, self.coef_cov = fit.est[0, :1] / Nbar, fit.jknife_se[0, :1] / Nbar, fit.jknife_cov[:1, :1] / Nbar**2
            self.intercept, self.intercept_se = float(fit.est[0, -1]), float(fit.jknife_se[0, -1])
            self.tot_delete_values = fit.delete_values[:, :1] @ M.T / Nbar
            self.intercept_delete_values = fit.delete_values[:, -1]

        def _combine_twostep_jknives(self, first, second, M_tot, c, Nbar=1):
            point = np.column_stack((second.est, first.est[:, -1]))
            deletes = np.column_stack((second.delete_values - c * (first.delete_values[:, -1] - first.est[0, -1]).reshape(-1, 1), first.delete_values[:, -1]))
            return SimpleNamespace(**jknife(point, deletes))

    class Hsq(MockBase):
        def _update_weights(self, x, w, N, M, tot, intercept):
            return 1 / (1 + x**2)

    class Gencov(Hsq):
        pass

    jk.LstsqJackknifeFast, jk.RatioJackknife = MockFast, MockRatio
    reg.IRWLS, reg.LD_Score_Regression, reg.Hsq, reg.Gencov = MockIRWLS, MockBase, Hsq, Gencov

    def mock_rg(data, args, log, M, refs, weight_name, index):
        events.append({'stock_received_rows': len(data), 'stock_received_chisq_max': args.chisq_max})
        final = data if args.chisq_max is None else data[data.Z1**2 * data.Z2**2 < args.chisq_max**2]
        shape = lambda col: final[col].to_numpy().reshape(-1, 1)
        x, w, n1, n2, z1, z2 = shape('L2'), shape('W_L2'), shape('N1'), shape('N2'), shape('Z1'), shape('Z2')
        masks = [None] * 3 if args.two_step is None else [z1**2 < 30, z2**2 < 30, (z1**2 < 30) & (z2**2 < 30)]
        first, second, cov = Hsq(z1**2, x, w, n1, M, 200, step1_ii=masks[0]), Hsq(z2**2, x, w, n2, M, 200, step1_ii=masks[1]), Gencov(z1*z2, x, w, np.sqrt(n1*n2), M, 200, step1_ii=masks[2])
        ratio = cov.tot / np.sqrt(first.tot * second.tot)
        fit = jk.RatioJackknife(np.array([[ratio]]), cov.tot_delete_values, np.sqrt(first.tot_delete_values * second.tot_delete_values))
        return SimpleNamespace(rg_ratio=ratio, rg_jknife=float(fit.jknife_est.item()), rg_se=float(fit.jknife_se.item()), _negative_hsq=None)

    ss = SimpleNamespace(_rg=mock_rg)
    originals = (ss._rg, MockBase.__init__, reg.IRWLS, MockFast.__init__, MockRatio.__init__, MockBase._combine_twostep_jknives, Hsq._update_weights, Gencov._update_weights)
    with tempfile.TemporaryDirectory(prefix='canonical_mock_fixture_') as temporary:
        root = Path(temporary)
        interval = root / 'intervals.tsv'
        with interval.open('w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['block_id', 'chromosome', 'start_inclusive', 'end_exclusive'])
            writer.writerows((b, 1, 10*b+1, 10*b+11) for b in range(200))
        blocks = {'build': 'GRCh37/hg19', 'block_count': 200, 'interval_tsv_path': str(interval), 'interval_tsv_sha256': sha(interval), 'reference_sha256': {str(i): 'a'*64 for i in range(1, 23)}}
        mock_plan = root / 'plan.json'
        mock_plan.write_text('{}')
        jobs = []
        for control in ('fixture_control_1', 'fixture_control_2'):
            for method in ('one_step', 'two_step30'):
                for partition in ('stock_default200', 'canonical200'):
                    for duplicate in ('A', 'B'):
                        identity = '__'.join((control, method, partition, duplicate))
                        out = root / identity
                        out.mkdir()
                        adapter = NativeCalibrationAdapter(blocks, partition, out / 'weighted_capture', np, admitted_scope='METHOD_CALIBRATION_ONLY')
                        adapter.coordinates = lambda wanted: {s: coordinates[s] for s in wanted}
                        args = SimpleNamespace(chisq_max=None, n_blocks=200, two_step=None if method == 'one_step' else 30)
                        with adapter.install(ss, reg, jk):
                            result = ss._rg(rows, args, None, np.array([[1000000.]]), ['L2'], 'W_L2', 0)
                        assert originals == (ss._rg, MockBase.__init__, reg.IRWLS, MockFast.__init__, MockRatio.__init__, MockBase._combine_twostep_jknives, Hsq._update_weights, Gencov._update_weights)
                        capture = adapter.record()
                        if partition == 'canonical200' and method == 'two_step30':
                            # Deliberately perturb a MOCK captured response by
                            # binary64-scale noise to exercise unequal-hash,
                            # numerically identical cross-arm comparisons.
                            # Both duplicates receive the same tiny fixture
                            # perturbation. No native captured data are edited.
                            item = capture['regressions'][0]['stages'][1]['final_weighted_response']
                            values = np.load(item['path'], allow_pickle=False)
                            values[0, 0] += 1e-14
                            with Path(item['path']).open('wb') as f:
                                np.save(f, values, allow_pickle=False)
                            item['sha256'] = sha(item['path'])
                        capture['point_estimates'] = adapter.attributes(result, ('rg_ratio', 'rg_jknife', 'rg_se', '_negative_hsq'))
                        job = {'job_id': identity, 'control_id': control, 'method': method, 'partition': partition, 'duplicate': duplicate, 'output_dir': str(out)}
                        jobs.append(job)
                        write_new(out / 'native_capture_receipt.json', {'job': job, 'plan_sha256': sha(mock_plan), 'capture': capture})
        plan = {'jobs': jobs, 'allowed_control_ids': ['fixture_control_1', 'fixture_control_2'], 'self_path': str(mock_plan)}
        arithmetic = validate(plan, np)
        assert arithmetic['all_checks_pass']
        unequal_hash_checks = [c for c in arithmetic['checks'] if c['name'].endswith('/hash_observation') and not c['binary_equal']]
        assert len(unequal_hash_checks) == 2
        # Force the last optional pair filter to remove exactly one row/block.
        # Stock still receives the unfiltered 600 rows and original argument.
        filtered_rows = rows.copy()
        filtered_rows.loc[filtered_rows.index % 3 == 0, ['Z1', 'Z2']] = 1000
        adapter = NativeCalibrationAdapter(blocks, 'canonical200', root / 'final_filter_capture', np, admitted_scope='METHOD_CALIBRATION_ONLY')
        adapter.coordinates = lambda wanted: {s: coordinates[s] for s in wanted}
        with adapter.install(ss, reg, jk):
            ss._rg(filtered_rows, SimpleNamespace(chisq_max=100, n_blocks=200, two_step=None), None, np.array([[1000000.]]), ['L2'], 'W_L2', 0)
        assert adapter.pair['ordered_final_snp_count'] == 400
        assert events[-1] == {'stock_received_rows': 600, 'stock_received_chisq_max': 100}
        # Fail an unsorted final set before any mock regression constructor and
        # verify all instrumented module bindings restore after the exception.
        adapter = NativeCalibrationAdapter(blocks, 'canonical200', root / 'bad_order_capture', np, admitted_scope='METHOD_CALIBRATION_ONLY')
        adapter.coordinates = lambda wanted: {s: coordinates[s] for s in wanted}
        try:
            with adapter.install(ss, reg, jk):
                ss._rg(rows.iloc[::-1], SimpleNamespace(chisq_max=None, n_blocks=200, two_step=None), None, np.array([[1000000.]]), ['L2'], 'W_L2', 0)
        except module.BoundaryAdmissionError as error:
            assert str(error) == 'FINAL_SNP_ORDER_NOT_GENOMIC'
        else:
            raise AssertionError('Bad order admitted')
        assert originals == (ss._rg, MockBase.__init__, reg.IRWLS, MockFast.__init__, MockRatio.__init__, MockBase._combine_twostep_jknives, Hsq._update_weights, Gencov._update_weights)
    receipt = {'schema': 'canonical_native_adapter_mock_fixture_v1', 'completed_utc': utc(), 'fixture_sha256': sha(__file__),
               'pinned_interface_checks': interfaces, 'mock_controls': 16, 'independent_arithmetic_checks': len(arithmetic['checks']),
               'all_mock_arithmetic_checks_pass': arithmetic['all_checks_pass'], 'original_final_filter_retained': True,
               'cross_arm_unequal_weighted_hash_numerical_fallback_checks': len(unequal_hash_checks),
               'mock_response_noise_only': 1e-14,
               'exception_restoration_pass': True, 'invented_SNP_count': 600, 'GWAS_inputs_read': False,
               'native_LDSC_modules_imported': False, 'native_estimator_calls': 0,
               'status': 'MOCK_INTERFACE_ARITHMETIC_PASS_ONLY_NOT_NATIVE_OR_REALISTIC_LD_VALIDATION'}
    write_new(PACKAGE / 'statistical_validation/canonical_native_calibration_fixture_receipt_v4_3.json', receipt)
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()

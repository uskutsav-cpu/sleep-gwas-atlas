"""Calibration-only scoped LDSC adapter. Pinned source files are never changed.

Stock _rg retains its own final filter and complete estimator control flow.
Only the canonical arm changes IRWLS jackknife separators. Instrumentation
records the exact final weighted design passed to the stock fast jackknife.
No native imports occur in this module; numerical modules are injected only
by the separately admitted worker. Not safe for concurrent fits/threads.
"""
from contextlib import contextmanager
import csv
import gzip
import hashlib
import importlib.util
from pathlib import Path

from canonical_calibration_common import PACKAGE, clean, sha


def prototype():
    path = PACKAGE / 'statistical_validation/common_boundary_adapter_prototype_v4.py'
    import sys
    if 'canonical_frozen_fixture_adapter' in sys.modules:
        module = sys.modules['canonical_frozen_fixture_adapter']
        if Path(module.__file__).resolve() != path.resolve():
            raise RuntimeError('FIXTURE_MODULE_BINDING_CHANGED')
        return module
    spec = importlib.util.spec_from_file_location('canonical_frozen_fixture_adapter', path)
    module = importlib.util.module_from_spec(spec)
    # dataclass consults sys.modules under some Python versions.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class NativeCalibrationAdapter:
    def __init__(self, block_receipt, partition, output, np, *, admitted_scope):
        if admitted_scope != 'METHOD_CALIBRATION_ONLY' or partition not in ('stock_default200', 'canonical200'):
            raise RuntimeError('NATIVE_ADAPTER_SCOPE_NOT_ADMITTED')
        self.receipt = block_receipt
        self.partition = partition
        self.output = Path(output)
        self.np = np
        self.output.mkdir(parents=True, exist_ok=False)
        module = prototype()
        module.verify_pinned_interfaces()
        self.module = module
        with Path(block_receipt['interval_tsv_path']).open(newline='') as f:
            self.intervals = [module.Interval(int(r['block_id']), int(r['chromosome']), int(r['start_inclusive']), int(r['end_exclusive'])) for r in csv.DictReader(f, delimiter='\t')]
        if sha(block_receipt['interval_tsv_path']) != block_receipt['interval_tsv_sha256'] or block_receipt['block_count'] != 200:
            raise RuntimeError('FROZEN_INTERVAL_IDENTITY_MISMATCH')
        self.adapter = None
        self.pair = None
        self.base = None
        self.stage = None
        self.regressions = []
        self.arrays = {}
        self.ratio = None

    def array(self, name, value, *, finite=True):
        a = self.np.asarray(value)
        if a.dtype.hasobject or (finite and not self.np.isfinite(a).all()):
            raise RuntimeError('INVALID_CAPTURE_ARRAY: ' + name)
        p = self.output / (name + '.npy')
        with p.open('xb') as f:
            self.np.save(f, a, allow_pickle=False)
        r = {'path': str(p), 'sha256': sha(p), 'shape': list(a.shape), 'dtype': str(a.dtype)}
        self.arrays[name] = r
        return r

    def attributes(self, obj, names):
        return {name: clean(getattr(obj, name, None), self.np) for name in names}

    def jackknife(self, obj):
        return self.attributes(obj, ('est', 'jknife_est', 'jknife_var', 'jknife_se', 'jknife_cov', 'delete_values', 'separators'))

    def coordinates(self, ids):
        wanted = set(ids)
        result = {}
        path = self.receipt['coordinate_map_path']
        if sha(path) != self.receipt['coordinate_map_sha256']:
            raise RuntimeError('COORDINATE_MAP_CHANGED')
        with gzip.open(path, 'rt', newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            if not {'SNP', 'CHR', 'BP'} <= set(reader.fieldnames):
                raise RuntimeError('COORDINATE_MAP_SCHEMA_CHANGED')
            for row in reader:
                if row['SNP'] in wanted:
                    if row['SNP'] in result:
                        raise RuntimeError('DUPLICATE_MAP_COORDINATE')
                    result[row['SNP']] = (int(row['CHR']), int(row['BP']))
        if set(result) != wanted:
            raise RuntimeError('FINAL_SNP_COORDINATE_UNAVAILABLE')
        return result

    @contextmanager
    def install(self, ss, reg, jk):
        np = self.np
        adapter = self
        old_rg = ss._rg
        old_base = reg.LD_Score_Regression.__init__
        old_irwls = reg.IRWLS
        old_fast = jk.LstsqJackknifeFast.__init__
        old_ratio = jk.RatioJackknife.__init__
        old_combine = reg.LD_Score_Regression._combine_twostep_jknives
        old_hsq_weights = reg.Hsq._update_weights
        old_cov_weights = reg.Gencov._update_weights

        def pair(rows, args, log, M, refs, weight_name, index):
            if adapter.pair is not None or len(adapter.regressions):
                raise RuntimeError('ONLY_ONE_FROZEN_PAIR_PER_WORKER')
            # Read-only identity calculation. Original _rg still executes its
            # exact pinned filtering branch using the ORIGINAL rows and args.
            final = rows if args.chisq_max is None else rows[rows.Z1**2 * rows.Z2**2 < args.chisq_max**2]
            ids = final.SNP.tolist()
            adapter.adapter = adapter.module.CanonicalAdapter(adapter.intervals, adapter.coordinates(ids), build=adapter.receipt['build'], reference_sha256=adapter.receipt['reference_sha256'])
            with adapter.adapter.pair_scope(ids, min(args.n_blocks, len(ids))):
                context = adapter.adapter.pair
                snp_path = adapter.output / 'final_ordered_SNP.txt.gz'
                with snp_path.open('xb') as raw:
                    with gzip.GzipFile(fileobj=raw, mode='wb', mtime=0) as out:
                        for snp in ids: out.write((snp + '\n').encode())
                adapter.pair = {'ordered_final_snp_count': len(ids), 'ordered_final_snp_sha256': hashlib.sha256(('\n'.join(ids) + '\n').encode()).hexdigest(),
                                'ordered_final_snp_file': str(snp_path), 'ordered_final_snp_file_sha256': sha(snp_path),
                                'block_counts': context['counts'], 'canonical_separators': context['separators'], 'args_two_step': args.two_step,
                                'final_filter_chisq_max': args.chisq_max, 'stock_final_filter_retained': True,
                                'N1': adapter.array('pair_N1', final.N1.to_numpy()), 'N2': adapter.array('pair_N2', final.N2.to_numpy()),
                                'raw_reference_design': adapter.array('pair_reference_LD', final[refs].to_numpy()),
                                'canonical_block_labels': adapter.array('pair_canonical_block_labels', np.array(context['labels'], dtype=np.int16))}
                result = old_rg(rows, args, log, M, refs, weight_name, index)
                if result is None or len(adapter.regressions) != 3 or adapter.ratio is None:
                    raise RuntimeError('NATIVE_CONTROL_FIT_FAILED_OR_RATIO_UNAVAILABLE')
                return result

        def base(obj, y, x, w, N, M, n_blocks, intercept=None, slow=False, step1_ii=None, old_weights=False):
            if adapter.pair is None or adapter.base is not None or n_blocks != 200 or slow or old_weights or intercept is not None or x.shape[1] != 1:
                raise RuntimeError('UNADMITTED_BASE_ESTIMATOR_PATH')
            names = ('hsq1', 'hsq2', 'gencov')
            name = names[len(adapter.regressions)]
            if type(obj).__name__ != ('Gencov' if name == 'gencov' else 'Hsq'):
                raise RuntimeError('REGRESSION_CALLBACK_ORDER_CHANGED')
            frame = {'name': name, 'type': type(obj).__name__, 'Nbar': float(np.mean(N)), 'M': clean(M, np), 'M_tot': float(np.sum(M)),
                     'stages': [], 'initial_weight_capture_count': 0, 'two_step': step1_ii is not None}
            if frame['two_step'] != (adapter.pair['args_two_step'] is not None):
                raise RuntimeError('METHOD_ARM_MASK_MISMATCH')
            if step1_ii is not None:
                frame['first_stage_row_indices'] = adapter.array(name + '_first_stage_row_indices', np.flatnonzero(step1_ii.reshape(-1)))
            adapter.base = frame
            try:
                with adapter.adapter.regression_scope(x.shape[0], step1_ii, old_weights):
                    old_base(obj, y, x, w, N, M, n_blocks, intercept=intercept, slow=slow, step1_ii=step1_ii, old_weights=old_weights)
                if frame['initial_weight_capture_count'] != 1:
                    raise RuntimeError('INITIAL_WEIGHT_CALLBACK_COUNT_CHANGED')
                frame['fitted_attributes'] = adapter.attributes(obj, ('tot', 'tot_se', 'tot_cov', 'coef', 'coef_se', 'coef_cov', 'intercept', 'intercept_se', 'twostep_filtered'))
                frame['combined_jackknife'] = adapter.jackknife(obj.jknife)
                frame['tot_delete_values'] = adapter.array(name + '_total_delete', obj.tot_delete_values)
                frame['intercept_delete_values'] = adapter.array(name + '_intercept_delete', obj.intercept_delete_values)
                if obj.twostep_filtered is not None and obj.twostep_filtered != x.shape[0] - frame['first_stage_row_indices']['shape'][0]:
                    raise RuntimeError('NATIVE_FIRST_STAGE_MASK_DISAGREES')
                adapter.regressions.append(frame)
            finally:
                adapter.base = None

        def weights(original):
            def captured(obj, *args, **kwargs):
                result = original(obj, *args, **kwargs)
                if adapter.base is None:
                    raise RuntimeError('INITIAL_WEIGHTS_OUTSIDE_BASE')
                adapter.base['initial_weight_capture_count'] += 1
                adapter.base['initial_weights'] = adapter.array(adapter.base['name'] + '_initial_weights', result)
                return result
            return captured

        class CapturedIRWLS(old_irwls):
            def __init__(obj, x, y, update_func, n_blocks, w=None, slow=False, separators=None):
                if adapter.base is None or adapter.stage is not None or n_blocks != 200 or slow:
                    raise RuntimeError('UNEXPECTED_IRWLS_CALLBACK')
                cursor = adapter.adapter.frame['cursor']
                prepared = adapter.adapter.frame['stages'][cursor]
                canonical = adapter.adapter.irwls_separators(x.shape[0], separators)
                stage = {'name': prepared[0], 'canonical_block_counts': prepared[3], 'canonical_separators': canonical,
                         'incoming_separators': clean(separators, np), 'incoming_separators_overridden': adapter.partition == 'canonical200', 'jackknife_calls': 0}
                adapter.stage = stage
                try:
                    super().__init__(x, y, update_func, n_blocks, w=w, slow=slow, separators=canonical if adapter.partition == 'canonical200' else separators)
                    if stage['jackknife_calls'] != 1:
                        raise RuntimeError('FINAL_WEIGHTED_JACKKNIFE_CALLBACK_COUNT_CHANGED')
                    adapter.base['stages'].append(stage)
                finally:
                    adapter.stage = None

        def fast(obj, x, y, n_blocks=None, separators=None):
            if adapter.stage is None or adapter.base is None:
                raise RuntimeError('FAST_JACKKNIFE_OUTSIDE_CAPTURED_IRWLS')
            prefix = adapter.base['name'] + '_' + adapter.stage['name']
            adapter.stage['jackknife_calls'] += 1
            adapter.stage['final_weighted_design'] = adapter.array(prefix + '_weighted_X', x)
            adapter.stage['final_weighted_response'] = adapter.array(prefix + '_weighted_y', y)
            old_fast(obj, x, y, n_blocks=n_blocks, separators=separators)
            cuts = np.asarray(obj.separators, dtype=np.int64)
            if len(cuts) != 201 or cuts[0] != 0 or cuts[-1] != x.shape[0] or not np.all(np.diff(cuts) > 0):
                raise RuntimeError('EMPTY_OR_INVALID_ACTUAL_JACKKNIFE_BLOCK')
            if adapter.partition == 'canonical200' and cuts.tolist() != adapter.stage['canonical_separators']:
                raise RuntimeError('CANONICAL_SEPARATORS_NOT_ACTUALLY_CONSUMED')
            adapter.stage['actual_separators'] = cuts.tolist()
            adapter.stage['actual_block_counts'] = np.diff(cuts).tolist()
            adapter.stage['jackknife'] = adapter.jackknife(obj)

        def combine(obj, first, second, M_tot, c, Nbar=1):
            result = old_combine(obj, first, second, M_tot, c, Nbar)
            if adapter.base is None or not adapter.base['two_step']:
                raise RuntimeError('TWO_STEP_COMBINATION_CONTEXT_MISMATCH')
            adapter.base['two_step_combination'] = {'c': clean(c, np), 'M_tot': clean(M_tot, np), 'Nbar': clean(Nbar, np),
                                                  'first_stage': adapter.jackknife(first), 'second_stage': adapter.jackknife(second), 'combined': adapter.jackknife(result)}
            return result

        def ratio(obj, est, numer_delete_values, denom_delete_values):
            old_ratio(obj, est, numer_delete_values, denom_delete_values)
            # Base._prop also uses RatioJackknife. Only capture the final RG.
            if adapter.pair is not None and adapter.base is None:
                if adapter.ratio is not None:
                    raise RuntimeError('UNEXPECTED_MULTIPLE_PAIR_RATIOS')
                adapter.ratio = {'native_jackknife': adapter.jackknife(obj), 'point_ratio': clean(est, np),
                                 'numerator_delete': adapter.array('ratio_numerator_delete', numer_delete_values),
                                 'denominator_delete': adapter.array('ratio_denominator_delete', denom_delete_values)}

        ss._rg = pair
        reg.LD_Score_Regression.__init__ = base
        reg.IRWLS = CapturedIRWLS
        jk.LstsqJackknifeFast.__init__ = fast
        jk.RatioJackknife.__init__ = ratio
        reg.LD_Score_Regression._combine_twostep_jknives = combine
        reg.Hsq._update_weights = weights(old_hsq_weights)
        reg.Gencov._update_weights = weights(old_cov_weights)
        try:
            yield
        finally:
            ss._rg = old_rg
            reg.LD_Score_Regression.__init__ = old_base
            reg.IRWLS = old_irwls
            jk.LstsqJackknifeFast.__init__ = old_fast
            jk.RatioJackknife.__init__ = old_ratio
            reg.LD_Score_Regression._combine_twostep_jknives = old_combine
            reg.Hsq._update_weights = old_hsq_weights
            reg.Gencov._update_weights = old_cov_weights

    def record(self):
        if self.pair is None or len(self.regressions) != 3 or self.ratio is None:
            raise RuntimeError('INCOMPLETE_NATIVE_CAPTURE')
        return {'partition': self.partition, 'pair': self.pair, 'regressions': self.regressions, 'ratio': self.ratio, 'arrays': self.arrays,
                'source_control_flow': 'Original _rg filter and estimator constructors retained; canonical arm replaces only final IRWLS separators.',
                'weighted_capture': 'Exact X and y passed from stock final IRWLS weighting to LstsqJackknifeFast; no weights or Nbar re-estimated by capture.',
                'biological_covariance_outcomes_computed': False, 'empirical_uncertainty_calibrated': False}

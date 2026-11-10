"""Pinned numerical utilities; avoid importing an executable trainer by basename."""
import numpy as np
from scipy.stats import spearmanr,pearsonr

def base_metrics(d, ids, prediction, medians, supported, pairs):
    cells = {}; yy = d['effect'][ids]
    for c in np.unique(d['cell'][ids]):
        mask = d['cell'][ids] == c; y = yy[mask]; p = prediction[mask]
        baseline = float(np.abs(y-medians[c]).mean()); error = float(np.abs(y-p).mean())
        cells[str(d['cell_names'][c])] = {'n': int(mask.sum()), 'MAE': error, 'training_cell_median_MAE': baseline, 'MAE_skill': 1-error/baseline if baseline else 0., 'magnitude_Spearman': float(spearmanr(y, p).statistic) if np.ptp(p) > 0 and np.ptp(y) > 0 else None, 'cell_supported_by_current_fit': bool(supported[c])}
    inverse = np.full(len(d['effect']), -1, dtype='int32'); inverse[ids] = np.arange(len(ids))
    comparisons = {}
    for label, pair in pairs.items():
        ok = np.all(inverse[pair] >= 0, axis=1); q = pair[ok]
        if not len(q):
            comparisons[label] = {'pairs': 0}; continue
        delta_y = d['effect'][q[:, 0]]-d['effect'][q[:, 1]]
        delta_p = prediction[inverse[q[:, 0]]]-prediction[inverse[q[:, 1]]]
        zero = float(np.abs(delta_y).mean()); mae = float(np.abs(delta_y-delta_p).mean())
        nonzero = delta_y != 0
        comparisons[label] = {'pairs': len(q), 'difference_MAE': mae, 'zero_difference_MAE': zero, 'skill_vs_zero_difference': 1-mae/zero if zero else 0., 'difference_Spearman': float(spearmanr(delta_y, delta_p).statistic) if np.ptp(delta_p) > 0 and np.ptp(delta_y) > 0 else None, 'ordering_accuracy_nonzero_measured_differences': float(np.mean(np.sign(delta_y[nonzero]) == np.sign(delta_p[nonzero]))) if nonzero.any() else None}
    same = [m['MAE_skill'] for m in cells.values() if m['cell_supported_by_current_fit']]
    return {'cells': cells, 'macro_MAE_skill_supported_cells': float(np.mean(same)) if same else None, 'comparisons': comparisons, 'all_raw_predictions_are_development_diagnostics': True}

def transformations(d, train, precision):
    nc = len(d['cell_names']); n = len(d['effect'])
    medians = np.full(nc, np.median(d['effect'][train]), dtype='float32'); supported = np.zeros(nc, dtype=bool)
    for c in np.unique(d['cell'][train]):
        medians[c] = np.median(d['effect'][train&(d['cell'] == c)]); supported[c] = True
    weights = np.ones(n, dtype='float32'); centers = np.zeros(n, dtype='float32'); quality = np.ones(n, dtype='float32'); metadata = {}
    contexts = np.unique(d['context'][train])
    for c in contexts:
        allrows = d['context'] == c; tr = train&allrows
        weights[allrows] = min(train.sum()/(len(contexts)*tr.sum()), 10.)
        a = tr&np.isfinite(d['activity_reference'])&np.isfinite(d['activity_alternate'])
        center = float(np.median(np.r_[d['activity_reference'][a], d['activity_alternate'][a]])) if a.any() else 0.
        centers[allrows] = center
        known = tr&np.isfinite(d['measurement_standard_error'])&(d['measurement_standard_error'] > 0)
        info = {'training_rows': int(tr.sum()), 'activity_center': center, 'activity_training_rows': int(a.sum()), 'SE_training_rows': int(known.sum())}
        if precision and known.any():
            tau = float(np.median(d['measurement_standard_error'][known])); allknown = allrows&np.isfinite(d['measurement_standard_error'])&(d['measurement_standard_error'] > 0)
            raw = np.clip(tau*tau/(tau*tau+d['measurement_standard_error'][allknown]**2), .1, 1.)
            trraw = np.clip(tau*tau/(tau*tau+d['measurement_standard_error'][known]**2), .1, 1.)
            normalization = float(trraw.mean()); quality[allknown] = np.clip(raw/normalization, .25, 2.)
            info.update(SE_scale=tau, precision_normalization=normalization)
        metadata[str(d['context_names'][c])] = info
    weights /= weights[train].mean()
    return medians, supported, weights, centers, quality, metadata

def detailed_metrics(d, ids, output, medians, supported, pairs):
    m = base_metrics(d, ids, output['magnitude'], medians, supported, pairs)
    inverse = np.full(len(d['effect']), -1, dtype='int32'); inverse[ids] = np.arange(len(ids))
    q = pairs['cross'][np.all(inverse[pairs['cross']] >= 0, axis=1)]
    if len(q):
        y = d['effect'][q[:, 0]]-d['effect'][q[:, 1]]
        p = output['magnitude'][inverse[q[:, 0]]]-output['magnitude'][inverse[q[:, 1]]]
        b = medians[d['cell'][q[:, 0]]]-medians[d['cell'][q[:, 1]]]
        bmae = float(np.abs(y-b).mean())
        strata = d['study'][q[:, 0]].astype('int64')*10000+d['cell'][q[:, 0]]*100+d['cell'][q[:, 1]]
        rho = []
        for s in np.unique(strata):
            ix = strata == s
            if ix.sum() >= 30 and np.ptp(y[ix]) > 0 and np.ptp(p[ix]) > 0:
                rho.append(float(spearmanr(y[ix], p[ix]).statistic))
        m['comparisons']['cross'].update(training_cell_median_difference_MAE=bmae, skill_vs_training_cell_median_difference=1-float(np.abs(y-p).mean())/bmae, macro_within_fixed_study_cell_pair_Spearman=float(np.mean(rho)) if rho else None, fixed_cell_pair_strata=len(rho))
    m['activity_auxiliary'] = {}
    for k in ['activity_reference', 'activity_alternate']:
        valid = np.isfinite(d[k][ids]); x = d[k][ids][valid]; p = output[k][valid]
        m['activity_auxiliary'][k] = {'n': int(valid.sum()), 'Pearson': float(pearsonr(x, p).statistic) if len(x) > 10 and np.ptp(x) > 0 and np.ptp(p) > 0 else None, 'scope': 'Original measured allele activity; correlation is separate from allelic-effect or magnitude prediction.'}
    return m

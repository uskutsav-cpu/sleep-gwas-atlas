"""Finite activity/precision development fits; no fresh-study outcomes read."""
import argparse, time
from activity_effect_network import *
from train_development import metrics as base_metrics
from scipy.stats import spearmanr, pearsonr


def partitions(d, study, phase):
    train = (d['split'] == 0) if phase == 'development' else (d['split'] != 2)
    if study == 'none':
        test = d['split'] == (1 if phase == 'development' else 2)
    else:
        sid = list(d['study_names']).index(study)
        groups = np.unique(d['group'][d['study'] == sid])
        train &= (d['study'] != sid)&~np.isin(d['group'], groups)
        test = (d['study'] == sid)&((d['split'] != 2) if phase == 'development' else (d['split'] == 2))
    assert not set(d['group'][train])&set(d['group'][test])
    return train, test


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--recipe', choices=['activity', 'activity_precision', 'local_activity_precision'], required=True)
    parser.add_argument('--held-study', default='none'); parser.add_argument('--phase', choices=['development', 'final', 'confidence'], default='development')
    parser.add_argument('--seed', type=int, default=42); parser.add_argument('--epochs', type=int, default=8); parser.add_argument('--avi', action='store_true')
    args = parser.parse_args(); assert_previous_intact()
    if args.phase != 'development':
        lock = json.loads((OUT/'activity_round_final_settings_lock.json').read_text())
        assert args.recipe == lock['selected_recipe'] and args.epochs == lock['selected_epoch']
    if args.avi:
        assert (OUT/'activity_round_model_selection.json').exists(), 'Select the sequence recipe before the matched AVI experiment'
    plan = json.loads((OUT/'measured_activity_round_plan.json').read_text())
    source = DATA/'measured_training_arrays_with_activity_and_SE.npz'; assert sha(source) == plan['input_data_sha256']
    with np.load(source) as z: d = {k: z[k] for k in z.files}
    with np.load(DATA/'measured_pair_arrays.npz') as z: pairs = {k: z[k] for k in z.files}
    train_mask, test_mask = partitions(d, args.held_study, args.phase)
    tr, te = np.flatnonzero(train_mask), np.flatnonzero(test_mask)
    assert len(tr) > 1000 and len(te) > 100
    medians, supported, weights, centers, quality, context_metadata = transformations(d, train_mask, args.recipe != 'activity')
    h = np.load(DATA/'pretrained_pooled128.npy', mmap_mode='r')
    local = np.load(DATA/'pretrained_variant_local128.npy', mmap_mode='r') if args.recipe.startswith('local_') else None
    avi_bank = np.load(DATA/'avi_training_bank.npy') if args.avi else None; avi_metadata = None
    if args.avi:
        vals = avi_bank[np.unique(d['bank_index'][tr])]; known = vals[:, 2] == 1
        median = np.nanmedian(vals[known, :2], axis=0); filled = np.where(np.isfinite(vals[:, :2]), vals[:, :2], median)
        mean, sd = filled.mean(0), filled.std(0); sd[sd < 1e-5] = 1
        avi_metadata = {'median': median.tolist(), 'mean': mean.tolist(), 'sd': sd.tolist()}
    torch.manual_seed(args.seed); model = ActivityEffectNetwork(medians, args.recipe, args.seed, args.avi)
    optimizer = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad], lr=3e-4, weight_decay=.001)
    initial = {k: v.detach().clone() for k, v in model.state_dict().items()}
    suffix = 'AVI' if args.avi else 'sequence'
    name = f'round2_{args.recipe}_{args.held_study}_{args.phase}_{suffix}_{args.seed}'
    assert not (OUT/(name+'_record.json')).exists(), 'Completed fit exists'
    info = {'recipe': args.recipe, 'held_study': args.held_study, 'phase': args.phase, 'seed': args.seed, 'AVI': args.avi, 'epochs': args.epochs, 'training_rows': len(tr), 'test_rows': len(te), 'training_indices_sha256': hashlib.sha256(tr.tobytes()).hexdigest(), 'test_indices_sha256': hashlib.sha256(te.tobytes()).hexdigest(), 'data_sha256': sha(source), 'activity_context_transformations': context_metadata, 'training_cell_medians': medians.tolist(), 'supported_cell_groups': supported.tolist(), 'AVI_transform': avi_metadata, 'source_code_sha256': {p.name: sha(p) for p in [Path(__file__), WORK/'activity_effect_network.py']}, 'fresh_study_performance_read': False, 'PID': os.getpid()}
    write_json(OUT/(name+'_started.json'), info)
    tp = {k: p[np.all(train_mask[p], axis=1)] for k, p in pairs.items()}
    def forward(ids):
        banks = d['bank_index'][ids]; values = torch.from_numpy(np.array(h[banks], copy=True)); cc = torch.from_numpy(d['cell'][ids].astype('int64'))
        loc = torch.from_numpy(np.array(local[banks], copy=True)) if local is not None else None
        av = None
        if args.avi:
            a = np.array(avi_bank[banks], copy=True); a[:, :2] = (np.where(np.isfinite(a[:, :2]), a[:, :2], median)-mean)/sd
            av = torch.from_numpy(np.clip(a, -12, 12))
        return model(values[:, 0], values[:, 1], cc, loc[:, 0] if loc is not None else None, loc[:, 1] if loc is not None else None, av)
    def predict():
        model.eval(); output = {k: np.empty(len(te), dtype='float32') for k in ['magnitude', 'signed_effect', 'activity_reference', 'activity_alternate']}
        with torch.no_grad():
            for i in range(0, len(te), 2048):
                got = forward(te[i:i+2048])
                for k in output: output[k][i:i+2048] = got[k].numpy()
        return output
    rng = np.random.default_rng(args.seed); best = -np.inf; curve = []; updates = 0; start = time.monotonic()
    checkpoint = MODELS/(name+'.safetensors')
    for epoch in range(1, args.epochs+1):
        torch.manual_seed(args.seed+epoch); model.train(); order = rng.permutation(tr); losses = []
        for i in range(0, len(order), 1024):
            main_ids = order[i:i+1024]; chunks = [main_ids]; counts = []
            for key in ['local', 'cross']:
                p = tp[key]; count = 32 if len(p) else 0; counts.append(count)
                if count: chunks.append(p[rng.integers(len(p), size=count)].reshape(-1))
            ids = np.concatenate(chunks); optimizer.zero_grad(); output = forward(ids)
            batch = {k: torch.from_numpy(d[k][ids]) for k in ['effect', 'signed_effect', 'activity_reference', 'activity_alternate']}
            batch['activity_center'] = torch.from_numpy(centers[ids])
            loss = activity_loss(output, batch, torch.from_numpy(weights[ids]), torch.from_numpy(quality[ids]), len(main_ids), *counts)
            assert torch.isfinite(loss); loss.backward()
            torch.nn.utils.clip_grad_norm_([p for p in model.parameters() if p.requires_grad], 1.)
            optimizer.step(); updates += 1; losses.append(float(loss.detach()))
        row = {'epoch': epoch, 'training_objective': float(np.mean(losses)), 'elapsed_seconds': time.monotonic()-start}
        if args.phase == 'development':
            output = predict(); met = detailed_metrics(d, te, output, medians, supported, pairs)
            row.update(met); score = met['macro_MAE_skill_supported_cells']; score = -1e10 if score is None else score
            if score > best:
                best = score; save_file({k: v.detach().contiguous() for k, v in model.state_dict().items()}, checkpoint)
                write_json(OUT/(name+'_best_development.json'), {'epoch': epoch, 'metrics': met, 'checkpoint_sha256': sha(checkpoint)})
                np.savez_compressed(OUT/(name+'_predictions.npz'), observation_index=te, **output)
        else:
            score = None
        curve.append(row); write_json(OUT/(name+'_curve.json'), curve)
        print('MEASURED ACTIVITY FIT', name, epoch, 'size_error_skill', score, 'seconds', round(time.monotonic()-start, 1), flush=True)
    if args.phase != 'development':
        # Fixed epochs only. Calibration outcomes never choose a checkpoint.
        save_file({k: v.detach().contiguous() for k, v in model.state_dict().items()}, checkpoint)
        output = predict(); np.savez_compressed(OUT/(name+'_predictions.npz'), observation_index=te, **output)
    saved = load_file(checkpoint)
    changes = {k: int(torch.count_nonzero(v != initial[k])) for k, v in saved.items()}
    write_json(OUT/(name+'_record.json'), {**info, 'actual_backprop_updates': updates, 'completed_epochs': args.epochs, 'checkpoint': str(checkpoint), 'checkpoint_sha256': sha(checkpoint), 'changed_scalar_counts': changes, 'magnitude_loss_precision_weighted': False, 'calibration_outcome_selected_epoch': False, 'scope': 'Development fit or fixed-epoch final/calibration fit as specified. This record does not establish useful outside prediction, confidence or a causal sleep–frailty mechanism.'})
    assert_previous_intact(); print('MEASURED ACTIVITY RUN COMPLETE', name, flush=True)


if __name__ == '__main__':
    main()

"""Real measured-magnitude fits; fresh Matoba effects are never read here."""
import argparse, time, collections
from small_effect_network import *
from scipy.stats import spearmanr


def metrics(d, ids, prediction, medians, supported, pairs):
    cells = {}; yy = d['effect'][ids]
    for c in np.unique(d['cell'][ids]):
        mask = d['cell'][ids] == c; y = yy[mask]; p = prediction[mask]
        baseline = float(np.abs(y-medians[c]).mean()); error = float(np.abs(y-p).mean())
        cells[str(d['cell_names'][c])] = {'n': int(mask.sum()), 'MAE': error, 'training_cell_median_MAE': baseline, 'MAE_skill': 1-error/baseline if baseline else 0., 'magnitude_Spearman': float(spearmanr(y, p).statistic) if np.std(p) > 1e-12 else 0., 'cell_supported_by_current_fit': bool(supported[c])}
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
        comparisons[label] = {'pairs': len(q), 'difference_MAE': mae, 'zero_difference_MAE': zero, 'skill_vs_zero_difference': 1-mae/zero if zero else 0., 'difference_Spearman': float(spearmanr(delta_y, delta_p).statistic) if np.std(delta_p) > 1e-12 else 0., 'ordering_accuracy_nonzero_measured_differences': float(np.mean(np.sign(delta_y[nonzero]) == np.sign(delta_p[nonzero]))) if nonzero.any() else None}
    same = [m['MAE_skill'] for m in cells.values() if m['cell_supported_by_current_fit']]
    return {'cells': cells, 'macro_MAE_skill_supported_cells': float(np.mean(same)) if same else None, 'comparisons': comparisons, 'all_raw_predictions_are_development_diagnostics': True}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--recipe', choices=['frozen', 'pooled_tune', 'end_to_end'], required=True)
    parser.add_argument('--held-study', default='none')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--epochs', type=int, default=6)
    parser.add_argument('--device', choices=['cpu', 'mps'], default='cpu')
    args = parser.parse_args()
    assert_previous_intact()
    fr = json.loads((OUT/'pretrained_feature_receipt.json').read_text())
    assert sha(DATA/'pretrained_pooled128.npy') == fr['sha256']
    with np.load(DATA/'measured_training_arrays.npz') as z:
        d = {k: z[k] for k in z.files}
    with np.load(DATA/'measured_pair_arrays.npz') as z:
        pairs = {k: z[k] for k in z.files}
    train_mask = d['split'] == 0
    if args.held_study == 'none':
        test_mask = d['split'] == 1
    else:
        sid = list(d['study_names']).index(args.held_study)
        excluded_groups = np.unique(d['group'][d['study'] == sid])
        train_mask &= (d['study'] != sid)&~np.isin(d['group'], excluded_groups)
        test_mask = (d['study'] == sid)&(d['split'] != 2)
        assert not set(d['group'][train_mask])&set(d['group'][test_mask])
    train = np.flatnonzero(train_mask); test = np.flatnonzero(test_mask)
    assert len(train) > 1000 and len(test) > 100
    assert not set(d['group'][train])&set(d['group'][test])
    nc = len(d['cell_names']); medians = np.full(nc, np.median(d['effect'][train]), dtype='float32')
    supported = np.zeros(nc, dtype=bool)
    for c in np.unique(d['cell'][train]):
        medians[c] = np.median(d['effect'][train[d['cell'][train] == c]])
        supported[c] = True
    weights = np.ones(len(d['effect']), dtype='float32')
    contexts = np.unique(d['context'][train])
    for c in contexts:
        ids = train[d['context'][train] == c]
        weights[d['context'] == c] = min(len(train)/(len(contexts)*len(ids)), 10.)
    weights /= weights[train].mean()
    train_pairs = {k: p[np.all(train_mask[p], axis=1)] for k, p in pairs.items()}
    h = np.load(DATA/'pretrained_pooled128.npy', mmap_mode='r')
    token_r = np.load(DATA/'tokens_reference150.npy', mmap_mode='r')
    token_a = np.load(DATA/'tokens_alternate150.npy', mmap_mode='r')
    torch.manual_seed(args.seed); np.random.seed(args.seed)
    model = SmallEffectNetwork(medians, args.recipe, args.seed).to(args.device)
    optimizer = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad], lr=3e-4, weight_decay=.001)
    initial = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
    name = f'{args.recipe}_{args.held_study}_{args.seed}'
    plan = {'recipe': args.recipe, 'held_study': args.held_study, 'seed': args.seed, 'fixed_epochs': args.epochs, 'training_rows': len(train), 'test_rows': len(test), 'training_indices_sha256': hashlib.sha256(train.tobytes()).hexdigest(), 'test_indices_sha256': hashlib.sha256(test.tobytes()).hexdigest(), 'source_features_sha256': fr['sha256'], 'data_sha256': sha(DATA/'measured_training_arrays.npz'), 'source_code_sha256': {'network': sha(WORK/'small_effect_network.py'), 'trainer': sha(Path(__file__))}, 'fresh_study_read': False, 'calibration_rows_used': False, 'held_source_neighborhoods_completely_removed_from_training': args.held_study != 'none', 'started_PID': os.getpid()}
    assert not (OUT/(name+'_record.json')).exists(), 'Completed fit exists; do not silently repeat'
    write_json(OUT/(name+'_started.json'), plan)
    if not (OUT/'full_model_search_started.json').exists():
        write_json(OUT/'full_model_search_started.json', {'initial_plan_sha256': sha(OUT/'initial_experiment_plan.json'), 'initial_recipes': ['frozen', 'pooled_tune', 'end_to_end'], 'first_run': name, 'fresh_study_performance_seen': False})
    def forward(ids):
        banks = d['bank_index'][ids]
        if args.recipe == 'end_to_end':
            unique, inverse = np.unique(banks, return_inverse=True)
            tokens = np.r_[token_r[unique], token_a[unique]].astype('int64')
            x = F.one_hot(torch.from_numpy(tokens), 4).permute(0, 2, 1).to(device=args.device, dtype=torch.float32)
            pooled = read_pooled(model.reader, x)
            ri = torch.from_numpy(inverse).to(args.device); ar, aa = pooled[:len(unique)][ri], pooled[len(unique):][ri]
        else:
            values = torch.from_numpy(np.array(h[banks], copy=True)).to(args.device)
            ar, aa = values[:, 0], values[:, 1]
        return model(ar, aa, torch.from_numpy(d['cell'][ids].astype('int64')).to(args.device))
    def predict():
        model.eval(); output = np.empty(len(test), dtype='float32')
        with torch.no_grad():
            for i in range(0, len(test), 2048 if args.recipe != 'end_to_end' else 128):
                n = 2048 if args.recipe != 'end_to_end' else 128
                output[i:i+n] = forward(test[i:i+n])['magnitude'].cpu().numpy()
        return output
    rng = np.random.default_rng(args.seed); curve = []; best = -np.inf; updates = 0; start = time.monotonic()
    baseline = metrics(d, test, medians[d['cell'][test]], medians, supported, pairs)
    write_json(OUT/(name+'_initial_median_comparison.json'), baseline)
    batch_size = 1024 if args.recipe != 'end_to_end' else 128
    pair_count = 32 if args.recipe != 'end_to_end' else 8
    for epoch in range(1, args.epochs+1):
        torch.manual_seed(args.seed+epoch); model.train(); shuffled = rng.permutation(train); losses = []
        for offset in range(0, len(shuffled), batch_size):
            main_ids = shuffled[offset:offset+batch_size]; chunks = [main_ids]; counts = []
            for key in ['local', 'cross']:
                p = train_pairs[key]; count = pair_count if len(p) else 0; counts.append(count)
                if count:
                    chunks.append(p[rng.integers(len(p), size=count)].reshape(-1))
            ids = np.concatenate(chunks); optimizer.zero_grad()
            prediction = forward(ids)
            batch = {k: torch.from_numpy(d[k][ids]).to(args.device) for k in ['effect', 'signed_effect', 'activity_reference', 'activity_alternate']}
            loss = measured_loss(prediction, batch, torch.from_numpy(weights[ids]).to(args.device), len(main_ids), *counts)
            assert torch.isfinite(loss)
            loss.backward()
            torch.nn.utils.clip_grad_norm_([p for p in model.parameters() if p.requires_grad], 1.)
            optimizer.step(); updates += 1; losses.append(float(loss.detach().cpu()))
        output = predict(); result = metrics(d, test, output, medians, supported, pairs)
        score = result['macro_MAE_skill_supported_cells']; score = -1e10 if score is None else score
        row = {'epoch': epoch, 'mean_training_objective': float(np.mean(losses)), 'elapsed_seconds': time.monotonic()-start, **result}
        curve.append(row); write_json(OUT/(name+'_curve.json'), curve)
        if score > best:
            best = score; state = {k: v.detach().cpu().contiguous() for k, v in model.state_dict().items()}
            checkpoint = MODELS/(name+'.safetensors'); save_file(state, checkpoint)
            changes = {k: int(torch.count_nonzero(v != initial[k])) for k, v in state.items()}
            np.savez_compressed(OUT/(name+'_development_predictions.npz'), observation_index=test, magnitude=output, measured_magnitude=d['effect'][test], original_split=d['split'][test], cell=d['cell'][test], group=d['group'][test])
            write_json(OUT/(name+'_best_development.json'), {'epoch': epoch, 'metrics': result, 'checkpoint_sha256': sha(checkpoint), 'changed_scalar_counts': changes, 'not_fresh_validation': True})
        print('DEVELOPMENT FIT', name, 'epoch', epoch, 'macro_error_skill', round(score, 5), 'elapsed_seconds', round(time.monotonic()-start, 1), flush=True)
    write_json(OUT/(name+'_record.json'), {**plan, 'actual_backprop_updates': updates, 'completed_epochs': args.epochs, 'best_development_score': best, 'best': json.loads((OUT/(name+'_best_development.json')).read_text()), 'checkpoint': str(checkpoint), 'checkpoint_sha256': sha(checkpoint), 'fresh_Matoba_performance_read': False, 'scope': 'Development-selected research fit. Whole-study tests, confidence checks, AVI comparison and fresh locked test remain required.'})
    assert_previous_intact()
    print('DEVELOPMENT RUN COMPLETE', name, flush=True)


if __name__ == '__main__':
    main()

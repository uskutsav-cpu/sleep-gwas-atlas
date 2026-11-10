"""Fixed-epoch matched final/omitted-source fits after development selection."""
import argparse, time
from activity_effect_network import *
from train_activity_round import transformations


def shared_initial_hash(model, avi):
    h = hashlib.sha256()
    for k, v in sorted(model.state_dict().items()):
        if avi and k == 'mag_hidden.weight':
            v = v[:, :-3]
        a = v.detach().cpu().numpy().copy()
        h.update(k.encode()); h.update(str(a.shape).encode()); h.update(a.tobytes())
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--held-study', default='none'); parser.add_argument('--seed', type=int, required=True); parser.add_argument('--avi', action='store_true')
    args = parser.parse_args(); assert_previous_intact()
    lock = json.loads((OUT/'activity_round_final_settings_lock.json').read_text())
    family, recipe, epochs = lock['selected_family'], lock['selected_recipe'], lock['selected_epoch']
    assert args.seed in lock['seeds'] and (args.held_study == 'none' or args.held_study in lock['source_omission_studies'])
    source = DATA/('measured_training_arrays.npz' if family == 'initial' else 'measured_training_arrays_with_activity_and_SE.npz')
    with np.load(source) as z: d = {k: z[k] for k in z.files}
    if family == 'initial': d['measurement_standard_error'] = np.full(len(d['effect']), np.nan, dtype='float32')
    with np.load(DATA/'measured_pair_arrays.npz') as z: pairs = {k: z[k] for k in z.files}
    train_mask = d['split'] != 2
    if args.held_study == 'none':
        test_mask = d['split'] == 2
    else:
        sid = list(d['study_names']).index(args.held_study); groups = np.unique(d['group'][d['study'] == sid])
        train_mask &= (d['study'] != sid)&~np.isin(d['group'], groups)
        test_mask = d['study'] == sid
    tr, te = np.flatnonzero(train_mask), np.flatnonzero(test_mask)
    assert not set(d['group'][tr])&set(d['group'][te])
    med, support, weights, centers, quality, context_metadata = transformations(d, train_mask, family == 'activity' and recipe != 'activity')
    nc = len(med); source_count = np.zeros(nc, dtype='int16')
    for c in np.unique(d['cell'][tr]): source_count[c] = len(np.unique(d['study'][tr[d['cell'][tr] == c]]))
    h = np.load(DATA/'pretrained_pooled128.npy', mmap_mode='r')
    local = np.load(DATA/'pretrained_variant_local128.npy', mmap_mode='r') if recipe.startswith('local_') else None
    avi_bank = np.load(DATA/'avi_training_bank.npy') if args.avi else None; avi_transform = None
    if args.avi:
        values = avi_bank[np.unique(d['bank_index'][tr])]; known = values[:, 2] == 1
        median = np.nanmedian(values[known, :2], axis=0); filled = np.where(np.isfinite(values[:, :2]), values[:, :2], median)
        mean, sd = filled.mean(0), filled.std(0); sd[sd < 1e-5] = 1
        avi_transform = {'median': median.tolist(), 'mean': mean.tolist(), 'sd': sd.tolist()}
    torch.manual_seed(args.seed)
    model = SmallEffectNetwork(med, recipe, args.seed, args.avi) if family == 'initial' else ActivityEffectNetwork(med, recipe, args.seed, args.avi)
    arm = 'AVI' if args.avi else 'sequence'; name = f'fixed_{family}_{recipe}_{args.held_study}_{arm}_{args.seed}'
    assert not (OUT/(name+'_record.json')).exists()
    checkpoint = MODELS/(name+'.safetensors')
    initial_hash = shared_initial_hash(model, args.avi)
    metadata = {'family': family, 'recipe': recipe, 'epoch_count': epochs, 'seed': args.seed, 'held_study': args.held_study, 'arm': arm, 'training_rows': len(tr), 'prediction_rows': len(te), 'training_cell_medians': med.tolist(), 'training_cell_source_counts': source_count.tolist(), 'cell_names': d['cell_names'].tolist(), 'context_transformations': context_metadata, 'AVI_transform': avi_transform, 'shared_initial_parameter_sha256': initial_hash, 'training_indices_sha256': hashlib.sha256(tr.tobytes()).hexdigest(), 'predict_indices_sha256': hashlib.sha256(te.tobytes()).hexdigest(), 'source_data_sha256': sha(source), 'settings_lock_sha256': sha(OUT/'activity_round_final_settings_lock.json'), 'code_sha256': sha(Path(__file__)), 'PID': os.getpid(), 'fresh_performance_used': False}
    write_json(OUT/(name+'_started.json'), metadata)
    optimizer = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad], lr=3e-4, weight_decay=.001)
    tp = {k: p[np.all(train_mask[p], axis=1)] for k, p in pairs.items()}
    def forward(ids):
        banks = d['bank_index'][ids]; v = torch.from_numpy(np.array(h[banks], copy=True)); cc = torch.from_numpy(d['cell'][ids].astype('int64'))
        av = None
        if args.avi:
            a = np.array(avi_bank[banks], copy=True); a[:, :2] = (np.where(np.isfinite(a[:, :2]), a[:, :2], median)-mean)/sd
            av = torch.from_numpy(np.clip(a, -12, 12))
        if family == 'initial': return model(v[:, 0], v[:, 1], cc, av)
        q = torch.from_numpy(np.array(local[banks], copy=True)) if local is not None else None
        return model(v[:, 0], v[:, 1], cc, q[:, 0] if q is not None else None, q[:, 1] if q is not None else None, av)
    rng = np.random.default_rng(args.seed); updates = 0; curve = []; start = time.monotonic()
    for epoch in range(1, epochs+1):
        torch.manual_seed(args.seed+epoch); model.train(); order = rng.permutation(tr); losses = []
        for i in range(0, len(order), 1024):
            main_ids = order[i:i+1024]; chunks = [main_ids]; counts = []
            for key in ['local', 'cross']:
                p = tp[key]; count = 32 if len(p) else 0; counts.append(count)
                if count: chunks.append(p[rng.integers(len(p), size=count)].reshape(-1))
            ids = np.concatenate(chunks); optimizer.zero_grad(); prediction = forward(ids)
            batch = {k: torch.from_numpy(d[k][ids]) for k in ['effect', 'signed_effect', 'activity_reference', 'activity_alternate']}
            w = torch.from_numpy(weights[ids])
            if family == 'initial':
                loss = measured_loss(prediction, batch, w, len(main_ids), *counts)
            else:
                batch['activity_center'] = torch.from_numpy(centers[ids])
                loss = activity_loss(prediction, batch, w, torch.from_numpy(quality[ids]), len(main_ids), *counts)
            assert torch.isfinite(loss); loss.backward()
            torch.nn.utils.clip_grad_norm_([p for p in model.parameters() if p.requires_grad], 1.)
            optimizer.step(); updates += 1; losses.append(float(loss.detach()))
        curve.append({'epoch': epoch, 'training_objective': float(np.mean(losses)), 'elapsed_seconds': time.monotonic()-start})
        write_json(OUT/(name+'_curve.json'), curve)
        print('FIXED CANDIDATE FIT', name, epoch, round(time.monotonic()-start, 1), flush=True)
    model.eval(); output = {k: np.empty(len(te), dtype='float32') for k in ['magnitude', 'signed_effect', 'activity_reference', 'activity_alternate']}
    with torch.no_grad():
        for i in range(0, len(te), 2048):
            got = forward(te[i:i+2048])
            for k in output: output[k][i:i+2048] = got[k].numpy()
    save_file({k: v.detach().contiguous() for k, v in model.state_dict().items()}, checkpoint)
    predictions = OUT/(name+'_predictions.npz')
    np.savez_compressed(predictions, observation_index=te, original_split=d['split'][te], **output)
    write_json(OUT/(name+'_record.json'), {**metadata, 'actual_backprop_updates': updates, 'checkpoint': str(checkpoint), 'checkpoint_sha256': sha(checkpoint), 'predictions': str(predictions), 'predictions_sha256': sha(predictions), 'completed_epochs': epochs, 'no_calibration_or_test_checkpoint_selection': True, 'scope': 'Fixed matched source-omitted or final candidate. Calibration and fresh-study performance remain separate checks.'})
    assert_previous_intact(); print('FIXED CANDIDATE COMPLETE', name, flush=True)


if __name__ == '__main__':
    main()

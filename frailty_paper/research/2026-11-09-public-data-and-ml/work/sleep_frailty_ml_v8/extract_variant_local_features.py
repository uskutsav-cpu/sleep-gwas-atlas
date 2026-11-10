"""Preserve pretrained activations near the actual substituted base."""
import time
from small_effect_network import *


def local_pool(activations, positions, downsample=16):
    centers = torch.div(positions, downsample, rounding_mode='floor').clamp(0, activations.shape[-1]-1)
    offsets = torch.tensor([-1, 0, 1], device=activations.device)
    ids = centers[:, None]+offsets[None, :]
    valid = (ids >= 0)&(ids < activations.shape[-1])
    safe = ids.clamp(0, activations.shape[-1]-1)
    values = torch.gather(activations, 2, safe[:, None, :].expand(-1, activations.shape[1], -1))
    return (values*valid[:, None, :]).sum(-1)/valid.sum(-1)[:, None]


def main():
    assert_previous_intact()
    assert (OUT/'measured_activity_round_plan.json').exists()
    with np.load(DATA/'measured_training_arrays.npz') as z:
        used = np.unique(z['bank_index'])
    reference = np.load(DATA/'tokens_reference150.npy', mmap_mode='r')
    alternate = np.load(DATA/'tokens_alternate150.npy', mmap_mode='r')
    global_bank = np.load(DATA/'pretrained_pooled128.npy', mmap_mode='r')
    device = 'mps' if torch.backends.mps.is_available() else 'cpu'
    model = base_model().to(device)
    assert all(p == 2 for p in json.loads((SOURCES/'legnet_models_example.cfg').read_text())['pool_sizes'])
    target = DATA/'pretrained_variant_local128.npy'
    output = np.lib.format.open_memmap(target, mode='w+', dtype='float32', shape=(len(reference), 2, 128))
    output[:] = np.nan
    start = time.monotonic(); error = 0.; batch = 256
    for offset in range(0, len(used), batch):
        ids = used[offset:offset+batch]
        rr, aa = reference[ids], alternate[ids]
        changes = rr != aa; assert np.all(changes.sum(1) == 1)
        positions = changes.argmax(1)
        tokens = np.r_[rr, aa].astype('int64')
        x = F.one_hot(torch.from_numpy(tokens), 4).permute(0, 2, 1).to(device=device, dtype=torch.float32)
        pos = torch.from_numpy(np.r_[positions, positions]).to(device)
        with torch.no_grad():
            activations = model.main(model.stem(x))
            global_values = activations.mean(-1).cpu().numpy()
            local_values = local_pool(activations, pos).cpu().numpy()
        gold = np.r_[global_bank[ids, 0], global_bank[ids, 1]]
        error = max(error, float(np.max(abs(global_values-gold))))
        assert error < 1e-5
        output[ids, 0] = local_values[:len(ids)]; output[ids, 1] = local_values[len(ids):]
        if offset % (batch*64) == 0:
            output.flush(); print('ACTUAL VARIANT-LOCAL DNA FEATURES', offset, '/', len(used), flush=True)
    output.flush(); assert np.isfinite(output[used]).all()
    # A separate NumPy gather/mean on actual saved activations checks index
    # clipping and valid-neighbor denominators independently of torch.gather.
    acts = activations.cpu().numpy(); pp = pos.cpu().numpy(); expected = []
    for a, p in zip(acts, pp):
        center = min(int(p)//16, a.shape[-1]-1)
        expected.append(a[:, max(0, center-1):min(a.shape[-1], center+2)].mean(1))
    gather_error = float(np.max(abs(np.array(expected)-local_values)))
    assert gather_error < 1e-5
    write_json(OUT/'variant_local_feature_receipt.json', {
        'output': str(target), 'sha256': sha(target), 'used_physical_pairs': len(used),
        'shape': list(output.shape), 'native_150bp_actual_alleles': True,
        'pretrained_checkpoint_sha256': sha(MODELS/'author_WTC11_original.safetensors'),
        'global_feature_replay_max_error': error, 'independent_NumPy_local_gather_error': gather_error,
        'pooling_rule': 'Known actual SNP position divided by the four 2-fold pooling stages; average containing bin and valid immediate neighbors.',
        'fresh_test_used': False, 'elapsed_seconds': time.monotonic()-start,
        'scope': 'Additional positional molecular features for the declared development comparison. No label or phenotype information is used to construct them.'
    })
    print('VARIANT-LOCAL FEATURES COMPLETE', flush=True)


if __name__ == '__main__':
    main()

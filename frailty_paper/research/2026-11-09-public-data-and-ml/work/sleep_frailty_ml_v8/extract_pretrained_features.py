"""Exact pooled author features for fast, measured-data transfer experiments."""
import time
from small_effect_network import *


def main():
    assert_previous_intact()
    receipt = json.loads((OUT/'measured_training_data_audit.json').read_text())
    assert sha(DATA/'measured_training_arrays.npz') == receipt['output_arrays_sha256']
    with np.load(DATA/'measured_training_arrays.npz') as z:
        used = np.unique(z['bank_index'])
    reference = np.load(DATA/'tokens_reference150.npy', mmap_mode='r')
    alternate = np.load(DATA/'tokens_alternate150.npy', mmap_mode='r')
    model = base_model()
    x = F.one_hot(torch.from_numpy(np.asarray(reference[used[:8]], dtype='int64')), 4).permute(0, 2, 1).float()
    with torch.no_grad():
        gold = model(x)
        pooled = read_pooled(model, x)
        replay = model.head[3](pooled_to_embedding(model, pooled)).squeeze(1)
    error = float((gold-replay).abs().max())
    assert error < 1e-5, error
    device = 'mps' if torch.backends.mps.is_available() else 'cpu'
    model.to(device)
    target = DATA/'pretrained_pooled128.npy'
    output = np.lib.format.open_memmap(target, mode='w+', dtype='float32', shape=(len(reference), 2, 128))
    output[:] = np.nan
    start = time.monotonic(); batch = 256
    for offset in range(0, len(used), batch):
        ids = used[offset:offset+batch]
        tokens = np.r_[reference[ids], alternate[ids]].astype('int64')
        x = F.one_hot(torch.from_numpy(tokens), 4).permute(0, 2, 1).to(device=device, dtype=torch.float32)
        with torch.no_grad():
            h = read_pooled(model, x).cpu().numpy()
        output[ids, 0] = h[:len(ids)]; output[ids, 1] = h[len(ids):]
        if offset % (batch*32) == 0:
            output.flush()
            print('PRETRAINED DNA FEATURES', offset, '/', len(used), 'elapsed_seconds', round(time.monotonic()-start, 1), flush=True)
    output.flush()
    assert np.isfinite(output[used]).all()
    write_json(OUT/'pretrained_feature_receipt.json', {
        'file': str(target), 'sha256': sha(target), 'shape': list(output.shape),
        'used_physical_pairs': len(used), 'unrepresented_rows_are_NaN': True,
        'activity_forward_vs_commuted_pooled_forward_max_error': error,
        'actual_150bp_native_DNA_not_fabricated_padding': True,
        'author_pretrained_checkpoint_sha256': sha(MODELS/'author_WTC11_original.safetensors'),
        'training_data_audit_sha256': sha(OUT/'measured_training_data_audit.json'),
        'device': device, 'elapsed_seconds': time.monotonic()-start,
        'fresh_test_DNA_or_outcomes_used': False,
        'scope': 'Frozen experimental DNA-activity features; usefulness remains to be judged on held studies. Mapper/head weights can subsequently be fine-tuned; full convolutional tuning requires live DNA forward passes.'
    })
    assert_previous_intact()
    print('PRETRAINED MEASURED-ACTIVITY FEATURES COMPLETE', flush=True)


if __name__ == '__main__':
    main()

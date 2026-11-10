"""Small CPU portability check for the existing model; no training or new test.

Only explicitly listed, hash-verified public DNA examples, open reader weights
and existing project head weights are loaded. This probe is not the full cloud
training pipeline and cannot continue a Mac GPU RNG stream on other hardware.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys

sys.dont_write_bytecode = True
import numpy as np
import torch
from safetensors.torch import load_file


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024**2), b''):
            digest.update(chunk)
    return digest.hexdigest()


def module_from_file(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def local_pool(a, position):
    centers = (position // 16).clamp(0, a.shape[-1] - 1)
    indices = centers[:, None] + torch.tensor([-1, 0, 1], device=a.device)
    valid = (indices >= 0) & (indices < a.shape[-1])
    indices = indices.clamp(0, a.shape[-1] - 1)
    values = torch.gather(a, 2, indices[:, None, :].expand(-1, a.shape[1], -1))
    return (values * valid[:, None, :]).sum(-1) / valid.sum(-1)[:, None]


def summarize(g, l, p, n):
    gr, ga = g[:n], g[n:]; lr, la = l[:n], l[n:]; pr, pa = p[:n], p[n:]
    gd, ld = ga - gr, la - lr
    def cosine(a, b):
        return 1 - (a * b).sum(1) / torch.sqrt(((a * a).sum(1) * (b * b).sum(1)).clamp_min(1e-24))
    scalars = torch.stack([pr, pa, pa - pr, abs(pa - pr), abs(gd).mean(1),
        torch.sqrt((gd * gd).mean(1).clamp_min(1e-12)), cosine(gr, ga), abs(ld).mean(1),
        torch.sqrt((ld * ld).mean(1).clamp_min(1e-12)), cosine(lr, la),
        (gr.mean(1) + ga.mean(1)) * .5, (gr.std(1, correction=0) + ga.std(1, correction=0)) * .5,
        (lr.mean(1) + la.mean(1)) * .5, (lr.std(1, correction=0) + la.std(1, correction=0)) * .5,
        abs(gd).amax(1), abs(ld).amax(1)], 1)
    return scalars, torch.stack([(gr + ga) * .5, gd, (lr + la) * .5, ld], 1)


def probe(folder):
    config = json.loads((folder / 'native_probe_configuration.json').read_text())
    for filename, digest in config['inputs_sha256'].items():
        assert sha(folder / filename) == digest, filename
    architecture = module_from_file('probe_author_legnet', folder / 'legnet_model.py')
    head_module = module_from_file('probe_existing_head', folder / 'cell_matched_effect_network.py')
    reader = architecture.LegNet(in_ch=4, **config['reader_hyperparameters']).double().eval()
    reader.load_state_dict(load_file(str(folder / 'WTC11_author_test1_val2.safetensors')), strict=True)
    head = head_module.CellMatchedEffectHead(config['cell_medians'], config['line_medians'],
                                             config['line_parent'], 236).double().eval()
    weights = load_file(str(folder / 'frozen_reader_none_42_epoch3.safetensors'))
    assert all(name.startswith('head.') for name in weights)
    head.load_state_dict({name.removeprefix('head.'): value for name, value in weights.items()}, strict=True)
    reader.requires_grad_(False); head.requires_grad_(False)
    keys = ['reference_tokens', 'length', 'position', 'alt', 'bank', 'physical', 'channels', 'cell', 'line']
    with np.load(folder / 'public_DNA_probe_inputs.npz', allow_pickle=False) as archive:
        data = {key: archive[key] for key in keys}
    n = len(data['length']); offset = np.r_[0, np.cumsum(data['length'].astype('int64'))]
    projection = torch.from_numpy(np.load(folder / 'fixed_channel_projections.npy')[2].copy()).double()
    def predict(swap=False):
        dynamic = torch.zeros((n, 16), dtype=torch.float64)
        channels = torch.zeros((n, 4, 16), dtype=torch.float64)
        for length in np.unique(data['length']):
            ids = np.flatnonzero(data['length'] == length)
            ref = np.stack([data['reference_tokens'][offset[i]:offset[i+1]] for i in ids])
            alt = ref.copy(); alt[np.arange(len(ids)), data['position'][ids]] = data['alt'][ids]
            alleles = np.r_[alt, ref] if swap else np.r_[ref, alt]
            tokens = torch.nn.functional.one_hot(torch.from_numpy(alleles.astype('int64')), 4).permute(0, 2, 1).double()
            activation = reader.main(reader.stem(tokens))
            position = torch.from_numpy(np.r_[data['position'][ids], data['position'][ids]].astype('int64'))
            raw_readout = reader.head(reader.mapper(activation).mean(-1)).squeeze(-1)
            scalars, blocks = summarize(activation.mean(-1), local_pool(activation, position), raw_readout, len(ids))
            dynamic[ids] = scalars
            channels[ids] = (blocks @ projection) * torch.tensor([1., 1000., 1., 1000.])[None, :, None]
        fixed = torch.from_numpy(data['bank'].copy()).double()
        if swap:
            fixed[:, :, [0, 1]] = fixed[:, :, [1, 0]]; fixed[:, :, 2] *= -1
        raw = torch.cat([fixed[:, :2, :], dynamic[:, None, :]], 1)
        scalar = torch.cat([torch.cat([((raw[:, r, 0] + raw[:, r, 1]) * .5)[:, None], raw[:, r, 3:16]], 1) for r in range(3)], 1)
        old = torch.from_numpy(data['channels'].astype('float32')).double()
        projected = torch.cat([old[:, :2], channels[:, None]], 1)
        even = torch.stack([projected[:, :, 0, :], abs(projected[:, :, 1, :]),
                            projected[:, :, 2, :], abs(projected[:, :, 3, :])], 2)
        features = torch.cat([scalar, torch.from_numpy(data['physical']).double(), even.reshape(n, 192)], 1)
        features = ((features - torch.tensor(config['normalization']['mean'], dtype=torch.float64)) /
                    torch.tensor(config['normalization']['sd'], dtype=torch.float64)).clamp(-12, 12)
        cell, line = data['cell'], data['line']
        cs, ls = np.asarray(config['supported_cells']), np.asarray(config['supported_lines'])
        cf = np.eye(len(cs))[cell] * cs[cell, None]
        lf = np.eye(len(ls))[line] * ls[line, None]
        return head(features, raw[:, :, 0], raw[:, :, 1], torch.from_numpy(cell.astype('int64')),
                    torch.from_numpy(cf), torch.from_numpy(lf))
    with torch.no_grad():
        forward, reverse = predict(), predict(True)
    errors = {name: float(abs(forward[name].numpy()[:8] - np.asarray(values)).max())
              for name, values in config['saved_MPS_reference_predictions_first_eight'].items()}
    assert max(errors.values()) < config['absolute_tolerance'], errors
    assert torch.allclose(forward['magnitude'], reverse['magnitude'], atol=1e-12, rtol=0)
    assert torch.allclose(forward['signed_effect'], -reverse['signed_effect'], atol=1e-12, rtol=0)
    result = {'passed_numerical_portability_probe': True, 'platform': platform.platform(),
              'torch_version': torch.__version__, 'numpy_version': np.__version__, 'device': 'cpu_float64',
              'actual_DNA_examples': n, 'actual_DNA_lengths': np.unique(data['length']).tolist(),
              'reference_prediction_rows': 8, 'maximum_reference_errors': errors,
              'absolute_tolerance': config['absolute_tolerance'],
              'allele_swap_magnitude_invariance_and_signed_antisymmetry': True,
              'training_executed': False, 'new_study_validation_performed': False,
              'Mac_GPU_checkpoint_exact_continuation_claimed': False,
              'scope': 'Same existing fitted model on small historical DNA examples. This is environment portability, not improved accuracy, biological validation or the full cloud training pipeline.'}
    (folder / 'portable_probe_result.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--input-dir', type=Path, default=Path.cwd())
    args = parser.parse_args(); print(json.dumps(probe(args.input_dir)), flush=True)

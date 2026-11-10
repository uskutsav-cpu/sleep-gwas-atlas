"""Actual model/gradient/weight-update audit using measured development DNA."""
import time, copy
from legnet_runtime import *


def main():
    assert_previous_intact()
    start = time.monotonic()
    source = {r['variant_GRCh38']: r for r in csv_rows(V6/'outputs/kosicki_reconstructed_continuous_effects.csv')}
    chosen = []
    for r in csv_rows(V7/'sequence_bank.csv.gz'):
        if r['context_feature_source'] == 'kosicki' and r['split'] == 'train':
            old = source[r['variant_GRCh38']]
            assert int(r['changed_position']) == int(old['changed_position'])
            a, b = crop_pair(r)
            chosen.append((r, a, b, float(old['mean_log2_ratio_difference'])))
            if len(chosen) == 32:
                break
    assert len(chosen) == 32
    model = base_model()
    original = {k: v.clone() for k, v in model.state_dict().items()}
    x = encode([r[1] for r in chosen]+[r[2] for r in chosen])
    target = torch.tensor([r[3] for r in chosen])
    with torch.no_grad():
        cpu = model(x)
    device = 'mps' if torch.backends.mps.is_available() else 'cpu'
    model.to(device)
    with torch.no_grad():
        actual = model(x.to(device)).cpu()
    error = float((actual-cpu).abs().max())
    assert error < .002, error
    model.to('cpu').double()
    for p in model.parameters():
        p.requires_grad_(False)
    weight = model.head[3].weight
    weight.requires_grad_(True)
    xd = x[:4].double()
    yd = target[:2].double()
    # Use one REF/ALT pair in each half, not four reference sequences.
    xd = torch.cat([x[:2], x[32:34]]).double()
    def objective():
        activity = model(xd)
        delta = activity[2:]-activity[:2]
        return ((delta-yd)**2).mean()
    loss = objective()
    loss.backward()
    index = int(weight.grad.abs().argmax())
    analytic = float(weight.grad.reshape(-1)[index])
    assert abs(analytic) > 1e-9
    value = float(weight.detach().reshape(-1)[index])
    eps = 1e-4
    with torch.no_grad():
        weight.reshape(-1)[index] = value+eps
        plus = float(objective())
        weight.reshape(-1)[index] = value-eps
        minus = float(objective())
        weight.reshape(-1)[index] = value
    finite = (plus-minus)/(2*eps)
    relative = abs(finite-analytic)/max(abs(finite), abs(analytic), 1e-12)
    assert relative < 1e-5, (finite, analytic, relative)
    model = base_model().to(device)
    # Batch-normalization statistics stay fixed in this engineering audit.
    # Actual learned weights in the last convolutional stage DO change.
    for name, p in model.named_parameters():
        p.requires_grad_(name.startswith(('main.3.', 'mapper.', 'head.')))
    model.eval()
    optimizer = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad], lr=1e-4, weight_decay=1e-3)
    xx = x.to(device); yy = target.to(device)
    losses = []
    for step in range(3):
        optimizer.zero_grad()
        activity = model(xx)
        delta = activity[32:]-activity[:32]
        loss = (delta.abs()-yy.abs()).abs().mean()+.25*(delta-yy).abs().mean()
        loss.backward()
        assert all(p.grad is None or torch.isfinite(p.grad).all() for p in model.parameters())
        torch.nn.utils.clip_grad_norm_([p for p in model.parameters() if p.requires_grad], 1.)
        optimizer.step()
        losses.append(float(loss.detach().cpu()))
    if device == 'mps':
        torch.mps.synchronize()
    state = {k: v.detach().cpu().contiguous() for k, v in model.state_dict().items()}
    changes = {k: int(torch.count_nonzero(v != original[k])) for k, v in state.items()}
    assert sum(v for k, v in changes.items() if k.startswith('main.3.')) > 0
    assert sum(v for k, v in changes.items() if k.startswith('stem.')) == 0
    checkpoint = MODELS/'engineering_smoke_only_not_selected.safetensors'
    save_file(state, checkpoint)
    receipt = {
        'passed': True, 'device': device, 'torch_version': torch.__version__,
        'Python_version': sys.version, 'source_architecture_sha256': sha(SOURCES/'legnet_model.py'),
        'pretrained_checkpoint_sha256': sha(MODELS/'author_WTC11_original.safetensors'),
        'total_parameter_scalars': sum(p.numel() for p in model.parameters()),
        'actual_trainable_parameter_scalars': sum(p.numel() for p in model.parameters() if p.requires_grad),
        'CPU_to_selected_device_forward_max_absolute_difference': error,
        'independent_finite_difference': {'analytic': analytic, 'finite_difference': finite, 'relative_error': relative},
        'measured_training_examples': len(chosen), 'source': 'Kosicki2025, V7 training partition only',
        'row_identity_sha256': hashlib.sha256('|'.join(r[0]['variant_GRCh38'] for r in chosen).encode()).hexdigest(),
        'actual_backprop_updates': 3, 'engineering_batch_losses': losses,
        'changed_scalar_count': sum(changes.values()),
        'changed_last_convolutional_stage_scalars': sum(v for k, v in changes.items() if k.startswith('main.3.')),
        'frozen_stem_changed_scalars': sum(v for k, v in changes.items() if k.startswith('stem.')),
        'smoke_checkpoint_sha256': sha(checkpoint), 'elapsed_seconds': time.monotonic()-start,
        'fresh_Matoba_outcome_performance_read': False,
        'scope': 'Engineering and measured-data gradient/weight-update check. Three updates on 32 training examples are not model development completion, useful held-study prediction or biological validation.'
    }
    write_json(OUT/'small_model_runtime_and_gradient_audit.json', receipt)
    assert_previous_intact()
    print('SMALL PRETRAINED MODEL AUDITED; ACTUAL CONVOLUTIONAL WEIGHTS UPDATED', json.dumps(receipt), flush=True)


if __name__ == '__main__':
    main()

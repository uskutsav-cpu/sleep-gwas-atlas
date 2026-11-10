"""Independent NumPy prediction/loss replay and fitted-gradient checks."""
from small_effect_network import *


def numpy_forward(state, reference, alternate, cell):
    def bn(x, prefix, eps=1e-5):
        return (x-state[prefix+'.running_mean'])/np.sqrt(state[prefix+'.running_var']+eps)*state[prefix+'.weight']+state[prefix+'.bias']
    def linear(x, prefix):
        return x@state[prefix+'.weight'].T+state[prefix+'.bias']
    def silu(x):
        return x/(1+np.exp(-x))
    def embed(x):
        x = bn(x, 'reader.mapper.block.0')
        x = x@state['reader.mapper.block.1.weight'][:, :, 0].T+state['reader.mapper.block.1.bias']
        return silu(bn(linear(x, 'reader.head.0'), 'reader.head.1'))
    a, b = embed(reference), embed(alternate)
    weight = state['cell_activity.weight'][cell]; bias = state['cell_activity.bias'][cell]
    ar = (a*weight).sum(1)+bias; aa = (b*weight).sum(1)+bias; signed = aa-ar
    mean = (a+b)*.5
    centered = mean-mean.mean(1, keepdims=True)
    normalized = centered/np.sqrt((centered**2).mean(1, keepdims=True)+1e-5)
    flags = np.eye(len(state['mag_bias']))[cell]
    features = np.c_[abs(b-a)*10., normalized, abs(signed)*10., flags]
    hidden = silu(linear(features, 'mag_hidden'))
    raw = state['mag_bias'][cell]+linear(hidden, 'mag_out')[:, 0]
    return {'magnitude': np.logaddexp(0, raw), 'signed_effect': signed, 'activity_reference': ar, 'activity_alternate': aa}


def numpy_loss(pred, batch, weights, nmain, nlocal, ncross):
    result = np.mean(weights[:nmain]*abs(pred['magnitude'][:nmain]-batch['effect'][:nmain]))
    signed = np.isfinite(batch['signed_effect']); signed[nmain:] = False
    if signed.any():
        result += .25*np.mean(weights[signed]*abs(pred['signed_effect'][signed]-batch['signed_effect'][signed]))
    active = np.isfinite(batch['activity_reference'])&np.isfinite(batch['activity_alternate']); active[nmain:] = False
    if active.any():
        delta = batch['activity_alternate'][active]-batch['activity_reference'][active]
        residual = pred['signed_effect'][active]-delta; error = abs(residual)
        huber = np.where(error < 1, .5*error**2, error-.5)
        result += .05*np.mean(weights[active]*huber)
    for count, offset in [(nlocal, nmain), (ncross, nmain+2*nlocal)]:
        if count:
            p = pred['magnitude'][offset:offset+2*count].reshape(count, 2)
            y = batch['effect'][offset:offset+2*count].reshape(count, 2)
            result += .2*np.mean(abs((p[:, 0]-p[:, 1])-(y[:, 0]-y[:, 1])))
    return float(result)


def main():
    name = 'pooled_tune_none_42'
    record = json.loads((OUT/(name+'_record.json')).read_text())
    assert sha(record['checkpoint']) == record['checkpoint_sha256']
    with np.load(DATA/'measured_training_arrays.npz') as z:
        d = {k: z[k] for k in z.files}
    with np.load(DATA/'measured_pair_arrays.npz') as z:
        pairs = {k: z[k] for k in z.files}
    train = d['split'] == 0
    signed = np.flatnonzero(train&np.isfinite(d['signed_effect']))[:4]
    active = np.flatnonzero(train&np.isfinite(d['activity_reference']))[:4]
    q = [p[np.all(train[p], axis=1)][:1].reshape(-1) for p in pairs.values()]
    ids = np.r_[signed, active, *q]
    assert len(ids) == 12
    h = np.load(DATA/'pretrained_pooled128.npy', mmap_mode='r')
    v = np.asarray(h[d['bank_index'][ids]], dtype='float64')
    model = SmallEffectNetwork(np.ones(len(d['cell_names']))*.1, 'pooled_tune').double().eval()
    model.load_state_dict(load_file(record['checkpoint']), strict=True)
    # load_state_dict casts the original fitted float32 weights exactly to
    # double; this audit neither modifies nor refits the saved checkpoint.
    state = {k: x.detach().numpy().astype('float64') for k, x in model.state_dict().items()}
    cell = d['cell'][ids].astype('int64')
    batch = {k: d[k][ids].astype('float64') for k in ['effect', 'signed_effect', 'activity_reference', 'activity_alternate']}
    weights = np.linspace(.5, 1.5, len(ids))
    x, y = torch.from_numpy(v[:, 0]), torch.from_numpy(v[:, 1]); c = torch.from_numpy(cell)
    actual = model(x, y, c)
    expected = numpy_forward(state, v[:, 0], v[:, 1], cell)
    errors = {k: float(np.max(abs(actual[k].detach().numpy()-expected[k]))) for k in expected}
    assert max(errors.values()) < 1e-8, errors
    tb = {k: torch.from_numpy(v) for k, v in batch.items()}; tw = torch.from_numpy(weights)
    loss = measured_loss(actual, tb, tw, 8, 1, 1)
    independent = numpy_loss(expected, batch, weights, 8, 1, 1)
    assert abs(float(loss.detach())-independent) < 1e-8
    loss.backward(); finite_checks = []
    named = dict(model.named_parameters())
    for key in ['reader.mapper.block.1.weight', 'mag_hidden.weight', 'mag_out.weight', 'cell_activity.weight', 'mag_bias']:
        p = named[key]
        assert p.grad is not None and torch.isfinite(p.grad).all()
        index = int(p.grad.abs().argmax()); analytic = float(p.grad.reshape(-1)[index])
        assert abs(analytic) > 1e-10
        original = float(p.detach().reshape(-1)[index]); eps = 1e-6
        with torch.no_grad():
            p.reshape(-1)[index] = original+eps
            plus = float(measured_loss(model(x, y, c), tb, tw, 8, 1, 1))
            p.reshape(-1)[index] = original-eps
            minus = float(measured_loss(model(x, y, c), tb, tw, 8, 1, 1))
            p.reshape(-1)[index] = original
        finite = (plus-minus)/(2*eps)
        relative = abs(finite-analytic)/max(abs(finite), abs(analytic), 1e-12)
        assert relative < 1e-4, (key, analytic, finite, relative)
        finite_checks.append({'parameter': key, 'analytic': analytic, 'finite_difference': finite, 'relative_error': relative})
    assert sha(record['checkpoint']) == record['checkpoint_sha256']
    assert_previous_intact()
    write_json(OUT/'fitted_effect_network_independent_audit.json', {
        'passed': True, 'actual_fitted_checkpoint': record['checkpoint'], 'checkpoint_sha256': record['checkpoint_sha256'],
        'NumPy_prediction_max_errors': errors, 'NumPy_vs_torch_loss_absolute_error': abs(float(loss.detach())-independent),
        'fitted_gradient_checks': finite_checks,
        'real_training_examples_include': ['Original signed Kosicki outcomes', 'Separate measured per-allele activity', 'Actual nearby-variant comparison', 'Actual matched-cell comparison'],
        'new_fitting_or_saved_weight_change': False,
        'scope': 'Independent arithmetic and fitted-gradient correctness, not prediction-performance or biological validation.'
    })
    print('FITTED PRETRAINED-EFFECT MODEL INDEPENDENTLY REPLAYED', flush=True)


if __name__ == '__main__':
    main()

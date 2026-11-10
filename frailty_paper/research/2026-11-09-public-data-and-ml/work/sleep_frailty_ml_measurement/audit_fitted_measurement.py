"""Independent numerical replay of the two fitted development heads.

This checks implementation, not predictive usefulness. No reserved outcomes read.
"""
from measurement_common import *
from measurement_inputs import load_data, normalized_inputs
from folded_effect_network import FoldedEffectHead, measurement_auxiliary, torch
from training_utilities import transformations
from activity_effect_network import activity_loss
from audit_first_fitted_network import numpy_loss
from safetensors.torch import load_file
from scipy.optimize import brentq
from scipy.special import ndtr, logsumexp


def numpy_forward(s, mean, difference, physical, reference, alternate, cell, cf, lf):
    def layer(x, name): return x @ s[name+'.weight'].T + s[name+'.bias']
    def silu(x): return x / (1 + np.exp(np.clip(-x, -600, 600)))
    flags = np.c_[physical, cf, lf]
    weights = s['cell_activity_weight'][cell] + lf @ s['line_activity_weight']
    bias = s['cell_activity_bias'][cell] + lf @ s['line_activity_bias']
    raw_mean = (reference + alternate) / 2
    activity = (raw_mean*weights).sum(1) + bias
    activity += layer(silu(layer(np.c_[mean, raw_mean, flags], 'mean_hidden')), 'mean_output')[:, 0]
    def branch(change):
        return layer(silu(layer(silu(layer(np.c_[mean, change, flags], 'delta_hidden')), 'delta_hidden2')), 'delta_output')[:, 0]
    delta = ((alternate-reference)*weights).sum(1) + (branch(difference)-branch(-difference))/2
    residual = layer(silu(layer(np.c_[mean, flags, activity], 'noise_hidden')), 'noise_output')[:, 0]
    sigma = np.maximum(np.logaddexp(0, s['cell_noise_bias'][cell] + lf@s['line_noise_offset'] + residual), 1e-4)
    magnitude = np.array([brentq(lambda t: ndtr((t-d)/v)-ndtr((-t-d)/v)-.5, 0, abs(d)+8*v, xtol=1e-14) for d, v in zip(delta, sigma)])
    return dict(magnitude=magnitude, signed_effect=delta, activity_reference=activity-delta/2,
                activity_alternate=activity+delta/2, assay_SD_component=sigma, latent_absolute_component=abs(delta))


def numpy_auxiliary(p, b, se, w, q, n):
    good = np.isfinite(se[:n]) & (se[:n] > 0)
    if not good.any(): return 0.
    sd = np.maximum(se[:n][good], .001)
    y, d = b['effect'][:n][good], p['signed_effect'][:n][good]
    z = np.c_[-.5*((y-d)/sd)**2, -.5*((y+d)/sd)**2]
    deviance = np.maximum(np.log(2)-logsumexp(z, axis=1), 0)
    consistency = np.sqrt(1+2*deviance)-1
    a = abs(np.log(p['assay_SD_component'][:n][good])-np.log(sd))
    huber = np.where(a < 1, .5*a*a, a-.5)
    weight = w[:n][good]*q[:n][good]
    return float(.05*np.mean(weight*consistency)+.05*np.mean(weight*huber))


def main():
    assert_prior_versions()
    d, pairs, bank, x, extra = load_data()
    train = d['split'] == 0
    _, _, w, centers, q, _ = transformations(d, train, True)
    signed = np.flatnonzero(train & np.isfinite(d['signed_effect']))[:4]
    rna = np.flatnonzero(train & np.isfinite(d['activity_reference']))[-4:]
    pair_blocks = [p[np.all(train[p], axis=1)][0] for p in [pairs['local'], pairs['cross']]]
    ids = np.r_[signed, rna, *pair_blocks]
    nmain = len(signed)+len(rna)
    batch = {k:d[k][ids].astype('float64') for k in ['effect', 'signed_effect', 'activity_reference', 'activity_alternate']}
    batch['activity_center'] = centers[ids].astype('float64')
    se = d['measurement_standard_error'][ids].astype('float64')
    weights, quality = w[ids].astype('float64'), q[ids].astype('float64')
    tb = {k:torch.from_numpy(v) for k,v in batch.items()}
    tw, tq, tse = map(torch.from_numpy, [weights, quality, se])
    results = []
    for recipe in ['observation_only', 'measurement_consistency']:
        stem = DEST/'development'/f'{recipe}_none_42'
        record = json.loads(Path(str(stem)+'_record.json').read_text())
        checkpoint = Path(str(stem)+'.safetensors')
        assert sha(checkpoint) == record['checkpoint_sha256']
        b = d['bank_index'][ids]
        mean, difference, physical = [v.astype('float64') for v in normalized_inputs(extra, x, b, record['normalization'])]
        raw = bank[b].astype('float64'); cell = d['cell'][ids].astype('int64'); line = d['line'][ids]
        cf = np.eye(len(record['cell_medians']))[cell]*np.array(record['supported_cells'])[cell,None]
        lf = np.eye(len(record['line_medians']))[line]*np.array(record['supported_lines'])[line,None]
        values = [mean, difference, physical, raw[:,:,0], raw[:,:,1], cell, cf, lf]
        tensors = [torch.from_numpy(v.copy()) for v in values]
        model = FoldedEffectHead(record['cell_medians'], record['line_medians'], record['line_parent']).double().eval()
        model.load_state_dict(load_file(str(checkpoint)), strict=True)
        state = {k:v.detach().numpy().copy() for k,v in model.state_dict().items()}
        expected = numpy_forward(state, *values)
        got = model(*tensors)
        errors = {k:float(np.max(abs(got[k].detach().numpy()-expected[k]))) for k in expected}
        assert max(errors.values()) < 1e-8, errors
        def total(prediction):
            loss = activity_loss(prediction, tb, tw, tq, nmain, 1, 1)
            if recipe == 'measurement_consistency': loss += measurement_auxiliary(prediction, tb, tse, tw, tq, nmain)
            return loss
        loss = total(got)
        independent = numpy_loss(expected, batch, weights, quality, nmain, 1, 1)
        if recipe == 'measurement_consistency': independent += numpy_auxiliary(expected, batch, se, weights, quality, nmain)
        loss_error = abs(float(loss.detach())-independent)
        assert loss_error < 1e-8, loss_error
        loss.backward(); checks = []
        # These parameters cannot change the deliberately detached activity
        # conditioner. Thus ordinary finite differences test the same derivative.
        for name in ['delta_hidden.weight','delta_hidden2.weight','delta_output.weight','noise_hidden.weight','noise_output.weight','cell_noise_bias','line_noise_offset']:
            parameter = dict(model.named_parameters())[name]
            assert parameter.grad is not None and torch.isfinite(parameter.grad).all()
            index = int(parameter.grad.abs().argmax()); analytic = float(parameter.grad.reshape(-1)[index])
            assert abs(analytic) > 1e-10, (recipe, name)
            original = float(parameter.detach().reshape(-1)[index]); epsilon = 2e-6
            alternatives = []
            for sign in [1, -1]:
                altered = {k:v.copy() for k,v in state.items()}
                altered[name].reshape(-1)[index] = original+sign*epsilon
                p = numpy_forward(altered, *values)
                objective = numpy_loss(p, batch, weights, quality, nmain, 1, 1)
                if recipe == 'measurement_consistency': objective += numpy_auxiliary(p, batch, se, weights, quality, nmain)
                alternatives.append(objective)
            finite = (alternatives[0]-alternatives[1])/(2*epsilon)
            relative = abs(finite-analytic)/max(abs(finite),abs(analytic),1e-12)
            assert relative < 1e-4, (recipe,name,analytic,finite,relative)
            checks.append(dict(parameter=name, analytic=analytic, independent_finite_difference=finite, relative_error=relative))
        swapped = list(tensors); swapped[1] = -tensors[1]; swapped[3], swapped[4] = tensors[4], tensors[3]
        reverse = model(*swapped)
        assert torch.allclose(got['signed_effect'], -reverse['signed_effect'], atol=1e-10, rtol=0)
        assert torch.allclose(got['magnitude'], reverse['magnitude'], atol=1e-10, rtol=0)
        identical = list(tensors); identical[1] = torch.zeros_like(tensors[1]); identical[4] = tensors[3]
        assert torch.equal(model(*identical)['signed_effect'], torch.zeros(len(ids),dtype=torch.float64))
        assert sha(checkpoint) == record['checkpoint_sha256']
        results.append(dict(recipe=recipe,checkpoint_sha256=record['checkpoint_sha256'],NumPy_forward_errors=errors,NumPy_loss_error=loss_error,fitted_gradient_checks=checks))
    assert_prior_versions()
    write_json(DEST/'fitted_measurement_independent_audit.json',dict(passed=True,models=results,
       checked_actual_measured_development_rows=len(ids),fitted_allele_swap_and_identical_DNA_invariants=True,
       reserved_outcomes_read=False,saved_weights_changed=False,
       detached_activity_conditioner_declared=True,gradient_scope='Signed-effect and noise branches; activity conditioner is intentionally detached.',
       scientific_predictive_validity_established=False))
    print('FITTED MEASUREMENT NUMERICAL AUDIT PASS',flush=True)


if __name__ == '__main__': main()

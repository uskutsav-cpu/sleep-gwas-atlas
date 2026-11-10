"""Shared DNA-activity reader with measured magnitude and paired-allele heads."""
from legnet_runtime import *
import torch.nn as nn
import torch.nn.functional as F


def read_pooled(model, tokens):
    return model.main(model.stem(tokens)).mean(-1)


def pooled_to_embedding(model, pooled):
    # With frozen BN running statistics, the author's mapper is affine.
    # Temporal mean and this affine map commute, preserving the model exactly
    # without caching the full per-position activation tensor.
    normalized = model.mapper.block[0](pooled)
    conv = model.mapper.block[1]
    mapped = F.linear(normalized, conv.weight.squeeze(-1), conv.bias)
    return model.head[:3](mapped)


class SmallEffectNetwork(nn.Module):
    def __init__(self, cell_medians, recipe, seed=42, avi=False):
        super().__init__()
        self.reader = base_model()
        self.recipe = recipe
        self.avi = avi
        self.nc = len(cell_medians)
        self.cell_activity = nn.Linear(256, self.nc)
        with torch.no_grad():
            self.cell_activity.weight.copy_(self.reader.head[3].weight.repeat(self.nc, 1))
            self.cell_activity.bias.copy_(self.reader.head[3].bias.repeat(self.nc))
        self.context_norm = nn.LayerNorm(256, elementwise_affine=False)
        common_dim = 513+self.nc
        self.mag_hidden = nn.Linear(common_dim+(3 if avi else 0), 128)
        self.mag_out = nn.Linear(128, 1)
        generator = torch.Generator().manual_seed(seed+3181)
        common_weights = torch.randn((128, common_dim), generator=generator)/np.sqrt(common_dim)
        with torch.no_grad():
            self.mag_hidden.weight[:, :common_dim].copy_(common_weights)
            if avi:
                self.mag_hidden.weight[:, common_dim:].zero_()
            self.mag_hidden.bias.zero_()
            self.mag_out.weight.zero_(); self.mag_out.bias.zero_()
        medians = torch.tensor(cell_medians, dtype=torch.float32).clamp_min(1e-4)
        self.mag_bias = nn.Parameter(torch.log(torch.expm1(medians)))
        for name, p in self.reader.named_parameters():
            allowed = recipe != 'frozen' and name.startswith(('mapper.', 'head.0.', 'head.1.'))
            if recipe == 'end_to_end' and name.startswith('main.3.'):
                allowed = True
            p.requires_grad_(allowed)

    def train(self, mode=True):
        super().train(mode)
        # Do not transfer source running means from a single mixed batch.
        self.reader.eval()
        return self

    def forward(self, pooled_reference, pooled_alternate, cell, avi_features=None):
        n = len(cell)
        z = pooled_to_embedding(self.reader, torch.cat([pooled_reference, pooled_alternate]))
        a, b = z[:n], z[n:]
        weights = self.cell_activity.weight[cell]
        ar = (a*weights).sum(1)+self.cell_activity.bias[cell]
        aa = (b*weights).sum(1)+self.cell_activity.bias[cell]
        signed = aa-ar
        features = [(b-a).abs()*10., self.context_norm((a+b)*.5), signed.abs().unsqueeze(1)*10., F.one_hot(cell, self.nc).float()]
        if self.avi:
            assert avi_features is not None
            features.append(avi_features)
        hidden = F.silu(self.mag_hidden(torch.cat(features, 1)))
        hidden = F.dropout(hidden, p=.1, training=self.training)
        magnitude = F.softplus(self.mag_bias[cell]+self.mag_out(hidden).squeeze(1))
        return {'magnitude': magnitude, 'signed_effect': signed, 'activity_reference': ar, 'activity_alternate': aa}


def measured_loss(prediction, batch, weights, nmain=None, nlocal=0, ncross=0):
    # Size prediction is the primary task. Source activity and known signed
    # effects are auxiliary real measurements, never teacher-model labels.
    nmain = len(weights) if nmain is None else nmain
    magnitude = (prediction['magnitude'][:nmain]-batch['effect'][:nmain]).abs()
    loss = (weights[:nmain]*magnitude).mean()
    signed = torch.isfinite(batch['signed_effect'])
    signed[nmain:] = False
    if signed.any():
        loss = loss+.25*(weights[signed]*(prediction['signed_effect'][signed]-batch['signed_effect'][signed]).abs()).mean()
    activity = torch.isfinite(batch['activity_reference'])&torch.isfinite(batch['activity_alternate'])
    activity[nmain:] = False
    if activity.any():
        # Pairwise activity difference cancels source-specific promoter offset.
        delta = batch['activity_alternate'][activity]-batch['activity_reference'][activity]
        loss = loss+.05*(weights[activity]*F.smooth_l1_loss(prediction['signed_effect'][activity], delta, reduction='none')).mean()
    for count, offset in [(nlocal, nmain), (ncross, nmain+2*nlocal)]:
        if count:
            p = prediction['magnitude'][offset:offset+2*count].reshape(count, 2)
            y = batch['effect'][offset:offset+2*count].reshape(count, 2)
            loss = loss+.2*((p[:, 0]-p[:, 1])-(y[:, 0]-y[:, 1])).abs().mean()
    return loss

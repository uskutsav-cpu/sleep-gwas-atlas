"""Measured two-allele supervision and optional actual variant-local features."""
from small_effect_network import *


class ActivityEffectNetwork(SmallEffectNetwork):
    def __init__(self, medians, recipe, seed=42, avi=False):
        assert recipe in ['activity', 'activity_precision', 'local_activity_precision']
        super().__init__(medians, 'pooled_tune', seed, avi=False)
        self.recipe = recipe
        self.use_local = recipe.startswith('local_')
        self.avi = avi
        common_dim = 513+self.nc
        extra = 2+(256 if self.use_local else 0)+(3 if avi else 0)
        old = self.mag_hidden
        self.mag_hidden = nn.Linear(common_dim+extra, 128)
        with torch.no_grad():
            self.mag_hidden.weight.zero_()
            self.mag_hidden.weight[:, :common_dim].copy_(old.weight)
            self.mag_hidden.bias.copy_(old.bias)

    def forward(self, global_reference, global_alternate, cell, local_reference=None, local_alternate=None, avi_features=None):
        n = len(cell)
        z = pooled_to_embedding(self.reader, torch.cat([global_reference, global_alternate]))
        a, b = z[:n], z[n:]
        weights = self.cell_activity.weight[cell]
        ar = (a*weights).sum(1)+self.cell_activity.bias[cell]
        aa = (b*weights).sum(1)+self.cell_activity.bias[cell]
        signed = aa-ar; mean = (a+b)*.5
        pieces = [(b-a).abs()*10., self.context_norm(mean), signed.abs()[:, None]*10., F.one_hot(cell, self.nc).to(dtype=z.dtype)]
        pieces += [((ar+aa)*.5)[:, None], torch.log1p(torch.sqrt((mean*mean).mean(1)+1e-12))[:, None]]
        if self.use_local:
            assert local_reference is not None and local_alternate is not None
            q = pooled_to_embedding(self.reader, torch.cat([local_reference, local_alternate]))
            pieces.append((q[n:]-q[:n]).abs()*10.)
        if self.avi:
            assert avi_features is not None
            pieces.append(avi_features)
        hidden = F.silu(self.mag_hidden(torch.cat(pieces, 1)))
        hidden = F.dropout(hidden, p=.1, training=self.training)
        magnitude = F.softplus(self.mag_bias[cell]+self.mag_out(hidden).squeeze(1))
        return {'magnitude': magnitude, 'signed_effect': signed, 'activity_reference': ar, 'activity_alternate': aa}


def activity_loss(prediction, batch, weights, quality, nmain, nlocal, ncross):
    # The main observable size-error target remains unmodified by SE.
    main_error = (prediction['magnitude'][:nmain]-batch['effect'][:nmain]).abs()
    result = (weights[:nmain]*main_error).mean()
    signed = torch.isfinite(batch['signed_effect']); signed[nmain:] = False
    if signed.any():
        result = result+.25*(weights[signed]*quality[signed]*(prediction['signed_effect'][signed]-batch['signed_effect'][signed]).abs()).mean()
    active = torch.isfinite(batch['activity_reference'])&torch.isfinite(batch['activity_alternate']); active[nmain:] = False
    if active.any():
        # Offsets are fitted from TRAINING examples of the original context.
        a = batch['activity_reference'][active]-batch['activity_center'][active]
        b = batch['activity_alternate'][active]-batch['activity_center'][active]
        ea = F.smooth_l1_loss(prediction['activity_reference'][active], a, reduction='none')
        eb = F.smooth_l1_loss(prediction['activity_alternate'][active], b, reduction='none')
        result = result+.05*(weights[active]*(ea+eb)*.5).mean()
        ed = F.smooth_l1_loss(prediction['signed_effect'][active], b-a, reduction='none')
        result = result+.05*(weights[active]*quality[active]*ed).mean()
    for count, offset in [(nlocal, nmain), (ncross, nmain+2*nlocal)]:
        if count:
            p = prediction['magnitude'][offset:offset+2*count].reshape(count, 2)
            y = batch['effect'][offset:offset+2*count].reshape(count, 2)
            q = quality[offset:offset+2*count].reshape(count, 2).prod(1).sqrt()
            result = result+.2*(q*((p[:, 0]-p[:, 1])-(y[:, 0]-y[:, 1])).abs()).mean()
    return result

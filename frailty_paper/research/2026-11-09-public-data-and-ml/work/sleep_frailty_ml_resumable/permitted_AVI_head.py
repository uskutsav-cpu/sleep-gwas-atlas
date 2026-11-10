"""Add score inputs without changing any existing head parameter or prediction."""
import copy

import torch
from torch import nn


def expand_head_for_scores(base_head, base_feature_count=236):
    """Insert three zero-weight columns before the head's existing context inputs.

    Both comparison arms use this same expansion. The scoreless arm supplies
    zero score values but the same availability feature. Construct on CPU before
    moving the whole model to its fitting device. This preserves CPU RNG state
    and every old parameter, including a previously fitted, nonconstant head.
    """
    old = base_head.hidden
    if old.weight.device.type != 'cpu':
        raise ValueError('Expand the common CPU initialization before moving to the fitting device.')
    if not 0 < base_feature_count < old.in_features:
        raise ValueError('The feature/context boundary must lie inside the input matrix.')
    expanded = copy.deepcopy(base_head)
    with torch.random.fork_rng(devices=[]):
        new = nn.Linear(old.in_features + 3, old.out_features,
                        bias=old.bias is not None, dtype=old.weight.dtype)
    with torch.no_grad():
        new.weight[:, :base_feature_count].copy_(old.weight[:, :base_feature_count])
        new.weight[:, base_feature_count:base_feature_count + 3].zero_()
        new.weight[:, base_feature_count + 3:].copy_(old.weight[:, base_feature_count:])
        if old.bias is not None:
            new.bias.copy_(old.bias)
    new.train(old.training)
    new.weight.requires_grad_(old.weight.requires_grad)
    if old.bias is not None:
        new.bias.requires_grad_(old.bias.requires_grad)
    expanded.hidden = new
    return expanded

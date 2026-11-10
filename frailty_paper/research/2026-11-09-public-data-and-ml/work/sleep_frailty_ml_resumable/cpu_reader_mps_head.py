"""Resource-only device partition of the unchanged native model.

The complete original forward and training-mode methods are inherited. Reader
and feature calculations stay on CPU; the existing head receives differentiable
device copies. The main training batch, model parameters, DNA lengths and loss
are unchanged. This class alone does not change a running experiment.
"""
from resumable_common import ROOT
from bounded_network import NativeAdaptationModel, torch


def transfer_head_inputs(module, inputs):
    destination = module.output.weight.device
    return tuple(value.to(destination) for value in inputs)


class CPUReaderMPSHead(NativeAdaptationModel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._head_device_bridge = self.head.register_forward_pre_hook(transfer_head_inputs)

    def to(self, device):
        target = torch.device(device)
        if target.type not in ['cpu', 'mps']:
            raise ValueError('This audited resource candidate supports CPU and the existing Mac head device only.')
        super().to('cpu')
        self.head.to(target)
        return self

    def anchor_loss(self):
        return super().anchor_loss().to(self.head.output.weight.device)

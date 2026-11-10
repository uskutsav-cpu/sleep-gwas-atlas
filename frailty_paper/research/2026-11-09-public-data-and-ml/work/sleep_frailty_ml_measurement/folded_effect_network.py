"""Couple signed allele effects to observed magnitude with explicit assay noise.

This is an experimental working observation model, not a confidence guarantee.
"""
import numpy as np
import torch
from torch import nn
import torch.nn.functional as F
torch.set_num_threads(2)
NORMAL_ABSOLUTE_MEDIAN=0.6744897501960817

class FoldedNormalMedian(torch.autograd.Function):
    @staticmethod
    def forward(ctx,delta,sigma):
        assert torch.all(sigma>0);low=torch.zeros_like(delta);high=delta.abs()+8*sigma
        for step in range(36):
            middle=(low+high)*.5;cdf=.5*(torch.erf((middle-delta)/(sigma*2**.5))+torch.erf((middle+delta)/(sigma*2**.5)))
            low=torch.where(cdf<.5,middle,low);high=torch.where(cdf>=.5,middle,high)
        median=(low+high)*.5;ctx.save_for_backward(delta,sigma,median);return median

    @staticmethod
    def backward(ctx,gradient):
        delta,sigma,median=ctx.saved_tensors;t=torch.tanh(median*delta/(sigma*sigma))
        return gradient*t,gradient*(median-delta*t)/sigma

def observed_magnitude_median(delta,sigma):return FoldedNormalMedian.apply(delta,sigma)

class FoldedEffectHead(nn.Module):
    def __init__(self,cell_medians,line_medians,line_parent,seed=42):
        super().__init__();nc,nl=len(cell_medians),len(line_medians);self.nc=nc;self.nl=nl
        self.cell_activity_weight=nn.Parameter(torch.tensor([[0.,0.,1.]]*nc));self.line_activity_weight=nn.Parameter(torch.zeros(nl,3));self.cell_activity_bias=nn.Parameter(torch.zeros(nc));self.line_activity_bias=nn.Parameter(torch.zeros(nl))
        self.delta_hidden=nn.Linear(96+96+2+nc+nl,64);self.delta_hidden2=nn.Linear(64,16);self.delta_output=nn.Linear(16,1)
        self.mean_hidden=nn.Linear(96+3+2+nc+nl,32);self.mean_output=nn.Linear(32,1)
        self.noise_hidden=nn.Linear(96+2+nc+nl+1,32);self.noise_output=nn.Linear(32,1)
        generator=torch.Generator().manual_seed(seed+26718)
        with torch.no_grad():
            for layer in [self.delta_hidden,self.delta_hidden2,self.mean_hidden,self.noise_hidden]:layer.weight.copy_(torch.randn(layer.weight.shape,generator=generator)/np.sqrt(layer.in_features));layer.bias.zero_()
            for layer in [self.delta_output,self.mean_output,self.noise_output]:layer.weight.zero_();layer.bias.zero_()
        cm=torch.tensor(cell_medians,dtype=torch.float32).clamp_min(1e-4)/NORMAL_ABSOLUTE_MEDIAN;lm=torch.tensor(line_medians,dtype=torch.float32).clamp_min(1e-4)/NORMAL_ABSOLUTE_MEDIAN;bias=torch.log(torch.expm1(cm));self.cell_noise_bias=nn.Parameter(bias);self.line_noise_offset=nn.Parameter(torch.log(torch.expm1(lm))-bias[torch.tensor(line_parent,dtype=torch.int64)])

    def forward(self,mean,difference,physical,raw_reference,raw_alternate,cell,cell_flags,line_flags):
        flags=torch.cat([physical,cell_flags,line_flags],1);weights=self.cell_activity_weight[cell]+line_flags@self.line_activity_weight;bias=self.cell_activity_bias[cell]+line_flags@self.line_activity_bias
        raw_mean=(raw_reference+raw_alternate)*.5;activity_mean=(raw_mean*weights).sum(1)+bias+self.mean_output(F.silu(self.mean_hidden(torch.cat([mean,raw_mean,flags],1)))).squeeze(1)
        def odd_branch(change):return self.delta_output(F.silu(self.delta_hidden2(F.silu(self.delta_hidden(torch.cat([mean,change,flags],1)))))).squeeze(1)
        signed=((raw_alternate-raw_reference)*weights).sum(1)+(odd_branch(difference)-odd_branch(-difference))*.5
        noise_residual=self.noise_output(F.silu(self.noise_hidden(torch.cat([mean,flags,activity_mean.detach()[:,None]],1)))).squeeze(1)
        sigma=F.softplus(self.cell_noise_bias[cell]+line_flags@self.line_noise_offset+noise_residual).clamp_min(1e-4)
        magnitude=observed_magnitude_median(signed,sigma)
        return {'magnitude':magnitude,'signed_effect':signed,'activity_reference':activity_mean-signed*.5,'activity_alternate':activity_mean+signed*.5,'assay_SD_component':sigma,'latent_absolute_component':signed.abs()}

def measurement_auxiliary(prediction,batch,known_SE,weights,quality,nmain):
    known=torch.isfinite(known_SE[:nmain])&(known_SE[:nmain]>0)
    if not known.any():return prediction['magnitude'].sum()*0
    sd=known_SE[:nmain][known].clamp_min(1e-3);delta=prediction['signed_effect'][:nmain][known];y=batch['effect'][:nmain][known];sigma=prediction['assay_SD_component'][:nmain][known]
    # A robust folded-Gaussian consistency penalty. Its square-root transform
    # bounds quadratic tail influence; it is not a calibrated likelihood claim.
    components=torch.stack([-.5*((y-delta)/sd)**2,-.5*((y+delta)/sd)**2],1);deviance=(np.log(2.)-torch.logsumexp(components,1)).clamp_min(0);consistency=torch.sqrt(1+2*deviance)-1
    precision=F.smooth_l1_loss(torch.log(sigma),torch.log(sd),reduction='none');w=weights[:nmain][known]*quality[:nmain][known]
    return .05*(w*consistency).mean()+.05*(w*precision).mean()

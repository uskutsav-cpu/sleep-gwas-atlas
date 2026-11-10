"""Small trainable molecular-effect head over three real, frozen cell readers."""
import numpy as np
import torch
from torch import nn
import torch.nn.functional as F

torch.set_num_threads(2)


class CellMatchedEffectHead(nn.Module):
    def __init__(self,cell_medians,line_medians,line_parent,nfeatures=66,seed=42):
        super().__init__();nc,nl=len(cell_medians),len(line_medians)
        self.cell_activity_weight=nn.Parameter(torch.tensor([[0.,0.,1.]]*nc));self.line_activity_weight=nn.Parameter(torch.zeros(nl,3))
        self.cell_activity_bias=nn.Parameter(torch.zeros(nc));self.line_activity_bias=nn.Parameter(torch.zeros(nl))
        self.hidden=nn.Linear(nfeatures+2+nc+nl,64);self.hidden2=nn.Linear(64,32);self.output=nn.Linear(32,1)
        generator=torch.Generator().manual_seed(seed+8711)
        with torch.no_grad():
            self.hidden.weight.copy_(torch.randn(self.hidden.weight.shape,generator=generator)/np.sqrt(self.hidden.in_features));self.hidden.bias.zero_()
            self.hidden2.weight.copy_(torch.randn(self.hidden2.weight.shape,generator=generator)/8.);self.hidden2.bias.zero_();self.output.weight.zero_();self.output.bias.zero_()
        cm=torch.tensor(cell_medians,dtype=torch.float32).clamp_min(1e-4);lm=torch.tensor(line_medians,dtype=torch.float32).clamp_min(1e-4)
        cb=torch.log(torch.expm1(cm));self.cell_bias=nn.Parameter(cb);self.line_offset=nn.Parameter(torch.log(torch.expm1(lm))-cb[torch.tensor(line_parent,dtype=torch.int64)])

    def forward(self,x,raw_reader_reference,raw_reader_alternate,cell,cell_flags,line_flags):
        weight=self.cell_activity_weight[cell]+line_flags@self.line_activity_weight;bias=self.cell_activity_bias[cell]+line_flags@self.line_activity_bias
        reference=(raw_reader_reference*weight).sum(1)+bias;alternate=(raw_reader_alternate*weight).sum(1)+bias;signed=alternate-reference
        inputs=torch.cat([x,signed.abs()[:,None]*10.,((reference+alternate)*.5)[:,None],cell_flags,line_flags],1)
        hidden=F.silu(self.hidden(inputs));hidden=F.dropout(hidden,p=.1,training=self.training);hidden=F.silu(self.hidden2(hidden))
        magnitude=F.softplus(self.cell_bias[cell]+line_flags@self.line_offset+self.output(hidden).squeeze(1))
        return {'magnitude':magnitude,'signed_effect':signed,'activity_reference':reference,'activity_alternate':alternate}

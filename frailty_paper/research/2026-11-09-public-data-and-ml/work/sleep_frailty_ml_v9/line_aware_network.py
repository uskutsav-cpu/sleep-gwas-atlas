"""Paired DNA effects conditional on explicit cell/culture identity."""
from v9_common import *
from legnet_runtime import base_model
from small_effect_network import pooled_to_embedding
import torch
import torch.nn as nn
import torch.nn.functional as F

torch.set_num_threads(2)


class LineAwareEffectNetwork(nn.Module):
    def __init__(self,cell_medians,line_medians,line_parent,recipe,seed=42,avi=False):
        super().__init__();assert recipe in ['frozen','tuned']
        self.reader=base_model();self.nc=len(cell_medians);self.nl=len(line_medians);self.avi=avi
        self.cell_activity=nn.Linear(256,self.nc);self.line_activity=nn.Linear(256,self.nl)
        with torch.no_grad():
            self.cell_activity.weight.copy_(self.reader.head[3].weight.repeat(self.nc,1));self.cell_activity.bias.copy_(self.reader.head[3].bias.repeat(self.nc))
            self.line_activity.weight.zero_();self.line_activity.bias.zero_()
        self.context_norm=nn.LayerNorm(256,elementwise_affine=False)
        dim=771+self.nc+self.nl
        self.mag_hidden=nn.Linear(dim+(3 if avi else 0),128);self.mag_out=nn.Linear(128,1)
        generator=torch.Generator().manual_seed(seed+3181)
        with torch.no_grad():
            self.mag_hidden.weight[:,:dim].copy_(torch.randn((128,dim),generator=generator)/np.sqrt(dim))
            if avi:self.mag_hidden.weight[:,dim:].zero_()
            self.mag_hidden.bias.zero_();self.mag_out.weight.zero_();self.mag_out.bias.zero_()
        cm=torch.tensor(cell_medians,dtype=torch.float32).clamp_min(1e-4);lm=torch.tensor(line_medians,dtype=torch.float32).clamp_min(1e-4)
        cb=torch.log(torch.expm1(cm));lb=torch.log(torch.expm1(lm))-cb[torch.tensor(line_parent,dtype=torch.int64)]
        self.mag_cell_bias=nn.Parameter(cb);self.mag_line_offset=nn.Parameter(lb)
        for name,p in self.reader.named_parameters():p.requires_grad_(recipe=='tuned' and name.startswith(('mapper.','head.0.','head.1.')))

    def train(self,mode=True):
        super().train(mode);self.reader.eval();return self

    def forward(self,global_reference,global_alternate,local_reference,local_alternate,cell,line,cell_flags,line_flags,avi_features=None):
        n=len(cell);z=pooled_to_embedding(self.reader,torch.cat([global_reference,global_alternate]));a,b=z[:n],z[n:]
        local=pooled_to_embedding(self.reader,torch.cat([local_reference,local_alternate]));la,lb=local[:n],local[n:]
        cw=self.cell_activity.weight[cell];cb=self.cell_activity.bias[cell]
        lw=line_flags@self.line_activity.weight;linebias=line_flags@self.line_activity.bias
        ar=(a*(cw+lw)).sum(1)+cb+linebias;aa=(b*(cw+lw)).sum(1)+cb+linebias;mean=(a+b)*.5
        features=[(b-a).abs()*10.,self.context_norm(mean),(aa-ar).abs()[:,None]*10.,((ar+aa)*.5)[:,None],torch.log1p(torch.sqrt((mean*mean).mean(1)+1e-12))[:,None],(lb-la).abs()*10.,cell_flags,line_flags]
        if self.avi:assert avi_features is not None;features.append(avi_features)
        hidden=F.silu(self.mag_hidden(torch.cat(features,1)));hidden=F.dropout(hidden,p=.1,training=self.training)
        bias=self.mag_cell_bias[cell]+line_flags@self.mag_line_offset
        magnitude=F.softplus(bias+self.mag_out(hidden).squeeze(1))
        return {'magnitude':magnitude,'signed_effect':aa-ar,'activity_reference':ar,'activity_alternate':aa}

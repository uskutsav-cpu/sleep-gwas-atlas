"""Only change in the AVI arm: three zero-initialized added input columns."""
from cell_matched_effect_network import CellMatchedEffectHead,torch,nn

class MatchedCoverageHead(CellMatchedEffectHead):
    def __init__(self,cell_medians,line_medians,line_parent,seed=42,avi=False):
        super().__init__(cell_medians,line_medians,line_parent,66,seed);self.avi=avi
        if avi:
            old=self.hidden;self.hidden=nn.Linear(old.in_features+3,old.out_features)
            with torch.no_grad():
                self.hidden.weight.zero_();self.hidden.weight[:,:66].copy_(old.weight[:,:66]);self.hidden.weight[:,69:].copy_(old.weight[:,66:]);self.hidden.bias.copy_(old.bias)

    def forward(self,x,raw_reader_reference,raw_reader_alternate,cell,cell_flags,line_flags,avi_features=None):
        if self.avi:
            assert avi_features is not None and avi_features.shape==(len(x),3);x=torch.cat([x,avi_features],1)
        else:assert avi_features is None
        return super().forward(x,raw_reader_reference,raw_reader_alternate,cell,cell_flags,line_flags)

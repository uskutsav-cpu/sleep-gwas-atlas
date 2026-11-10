"""Isolated, unselected score-input extension of the pinned native DNA model.

The running reader-adaptation comparison does not import this module. Its
forward computation is the pinned bounded implementation plus one post-scaling
concatenation. Only designated downloadable composite AVI values may enter.
This implementation alone does not select a model or authorize an outcome join.
"""
import copy

from bounded_network import NativeAdaptationModel, local_pool, summarize, torch, np
from permitted_AVI_features import transform
from permitted_AVI_head import expand_head_for_scores


class PermittedAVIModel(NativeAdaptationModel):
    def __init__(self, medians, line_medians, line_parent, norm, recipe,
                 score_scaling, use_scores, seed=42):
        super().__init__(medians, line_medians, line_parent, norm, recipe, seed)
        self.head = expand_head_for_scores(self.head)
        self.score_scaling = copy.deepcopy(score_scaling)
        self.use_scores = bool(use_scores)

    def forward(self,ids,d,packed,bank,x,extra,supported,supported_lines,avi_bank,swap=False):
        device=self.feature_mean.device;dtype=self.feature_mean.dtype;b=d['bank_index'][ids];unique=b;inverse=np.arange(len(b),dtype="int64")
        summaries=torch.zeros((len(unique),16),device=device,dtype=dtype);channels=torch.zeros((len(unique),4,16),device=device,dtype=dtype)
        for length in np.unique(packed.length[unique]):
            selected=np.flatnonzero(packed.length[unique]==length);which=unique[selected]
            all_summaries=[];all_channels=[]
            for offset in range(0,len(which),32):
                real=which[offset:offset+32];nreal=len(real)
                used=np.r_[real,np.repeat(real[-1],32-nreal)]
                tokens=packed.batch(used,device,dtype)
                if swap:tokens=torch.cat([tokens[32:],tokens[:32]])
                with torch.no_grad():prefix=self.reader.main[:3](self.reader.stem(tokens))
                activations=self.reader.main[3](prefix);g=activations.mean(-1)
                pos=torch.from_numpy(np.r_[packed.position[used],packed.position[used]].astype('int64')).to(device)
                local=local_pool(activations,pos);p=self.reader.head(self.reader.mapper(activations).mean(-1)).squeeze(-1)
                summary,blocks=summarize(g,local,p,32)
                all_summaries.append(summary[:nreal]);all_channels.append(((blocks@self.projection)*self.projection_multiplier)[:nreal])
            summary=torch.cat(all_summaries);projected=torch.cat(all_channels)
            index=torch.from_numpy(selected).to(device);summaries=summaries.index_copy(0,index,summary);channels=channels.index_copy(0,index,projected)
        take=torch.from_numpy(inverse).to(device);dynamic=summaries[take];projected=channels[take]
        fixed=torch.from_numpy(bank[b].copy()).to(device=device,dtype=dtype)
        if swap:
            fixed=fixed.clone();fixed[:,:,[0,1]]=fixed[:,:,[1,0]];fixed[:,:,2]=-fixed[:,:,2]
        raw=torch.cat([fixed[:,:2,:],dynamic[:,None,:]],1)
        scalar=torch.cat([torch.cat([((raw[:,r,0]+raw[:,r,1])*.5)[:,None],raw[:,r,3:16]],1) for r in range(3)],1)
        old=torch.from_numpy(extra[b].astype('float32')).to(device=device,dtype=dtype)
        all_projected=torch.cat([old[:,:2],projected[:,None]],1)
        even=torch.stack([all_projected[:,:,0,:],abs(all_projected[:,:,1,:]),all_projected[:,:,2,:],abs(all_projected[:,:,3,:])],2)
        physical=torch.from_numpy(x[b,48:50].copy()).to(device=device,dtype=dtype)
        features=torch.cat([scalar,physical,even.reshape(len(b),192)],1)
        features=((features-self.feature_mean)/self.feature_sd).clamp(-12,12)
        score_features=torch.from_numpy(transform(avi_bank[b],self.score_scaling,self.use_scores)).to(device=device,dtype=dtype)
        features=torch.cat([features,score_features],1)
        c=d['cell'][ids];line=d['line'][ids]
        cf=np.eye(len(supported),dtype='float32')[c]*supported[c,None]
        lf=np.eye(len(supported_lines),dtype='float32')[line]*supported_lines[line,None]
        return self.head(features,raw[:,:,0],raw[:,:,1],torch.from_numpy(c.astype('int64')).to(device),
            torch.from_numpy(cf).to(device=device,dtype=dtype),torch.from_numpy(lf).to(device=device,dtype=dtype))

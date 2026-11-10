"""Change actual DNA-reader weights while retaining matched, symmetric heads."""
import importlib.util
from bounded_common import *
import torch
from torch import nn
import torch.nn.functional as F
from safetensors.torch import load_file
from cell_matched_effect_network import CellMatchedEffectHead

torch.set_num_threads(2)


def author_reader():
    receipt=json.loads((PREP/'cell_reader_tensor_and_architecture_audit.json').read_text());path=V8/'sources/legnet_model.py'
    assert sha(path)==receipt['author_architecture_sha256']
    assert sha(PREP/'WTC11_author_test1_val2.safetensors')==receipt['checks']['WTC11']['tensor_file_sha256']
    spec=importlib.util.spec_from_file_location('native_adaptation_author_architecture',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    cfg=json.loads((PREP/'WTC11_author_config.json').read_text());keys=['stem_ch','stem_ks','ef_ks','ef_block_sizes','pool_sizes','resize_factor']
    model=module.LegNet(in_ch=4,**{k:cfg[k] for k in keys});model.load_state_dict(load_file(str(PREP/'WTC11_author_test1_val2.safetensors')),strict=True)
    return model.eval()


def local_pool(a,position):
    centers=(position//16).clamp(0,a.shape[-1]-1);indices=centers[:,None]+torch.tensor([-1,0,1],device=a.device)
    valid=(indices>=0)&(indices<a.shape[-1]);indices=indices.clamp(0,a.shape[-1]-1)
    values=torch.gather(a,2,indices[:,None,:].expand(-1,a.shape[1],-1))
    return (values*valid[:,None,:]).sum(-1)/valid.sum(-1)[:,None]


def summarize(g,l,p,n):
    gr,ga=g[:n],g[n:];lr,la=l[:n],l[n:];pr,pa=p[:n],p[n:];gd=ga-gr;ld=la-lr
    def cosine(a,b):return 1-(a*b).sum(1)/torch.sqrt(((a*a).sum(1)*(b*b).sum(1)).clamp_min(1e-24))
    scalars=torch.stack([pr,pa,pa-pr,abs(pa-pr),abs(gd).mean(1),torch.sqrt((gd*gd).mean(1).clamp_min(1e-12)),cosine(gr,ga),
       abs(ld).mean(1),torch.sqrt((ld*ld).mean(1).clamp_min(1e-12)),cosine(lr,la),(gr.mean(1)+ga.mean(1))*.5,
       (gr.std(1,correction=0)+ga.std(1,correction=0))*.5,(lr.mean(1)+la.mean(1))*.5,
       (lr.std(1,correction=0)+la.std(1,correction=0))*.5,abs(gd).amax(1),abs(ld).amax(1)],1)
    blocks=torch.stack([(gr+ga)*.5,gd,(lr+la)*.5,ld],1)
    return scalars,blocks


class NativeAdaptationModel(nn.Module):
    def __init__(self,medians,line_medians,line_parent,norm,recipe,seed=42):
        super().__init__();assert recipe in ['frozen_reader','adapt_last_stage'];self.recipe=recipe
        self.reader=author_reader();self.head=CellMatchedEffectHead(medians,line_medians,line_parent,236,seed)
        for name,p in self.reader.named_parameters():p.requires_grad_(recipe=='adapt_last_stage' and name.startswith('main.3.'))
        self.register_buffer('feature_mean',torch.tensor(norm['mean']));self.register_buffer('feature_sd',torch.tensor(norm['sd']))
        self.register_buffer('projection',torch.from_numpy(np.load(CHANNEL/'fixed_channel_projections.npy')[2].copy()))
        self.register_buffer('projection_multiplier',torch.tensor([1.,1000.,1.,1000.])[None,:,None])
        self.anchored={name:p.detach().clone() for name,p in self.reader.named_parameters() if name.startswith('main.3.')}

    def train(self,mode=True):
        super().train(mode);self.reader.eval();return self

    def anchor_loss(self):
        if self.recipe=='frozen_reader':return self.head.output.weight.sum()*0
        total=0.;denominator=0.
        for name,p in self.reader.named_parameters():
            if name in self.anchored:
                original=self.anchored[name].to(device=p.device,dtype=p.dtype);total=total+(p-original).square().sum();denominator=denominator+original.square().sum()
        return .01*total/denominator.clamp_min(1e-12)

    def forward(self,ids,d,packed,bank,x,extra,supported,supported_lines,swap=False):
        device=self.feature_mean.device;dtype=self.feature_mean.dtype;b=d['bank_index'][ids];unique=b;inverse=np.arange(len(b),dtype="int64")
        summaries=torch.zeros((len(unique),16),device=device,dtype=dtype);channels=torch.zeros((len(unique),4,16),device=device,dtype=dtype)
        for length in np.unique(packed.length[unique]):
            selected=np.flatnonzero(packed.length[unique]==length);which=unique[selected]
            all_summaries=[];all_channels=[]
            for offset in range(0,len(which),32):
                real=which[offset:offset+32];nreal=len(real)
                # Exact real-construct copies bound the batch dimension only.
                # No extension or padding of the biological sequence is used.
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
        c=d['cell'][ids];line=d['line'][ids]
        cf=np.eye(len(supported),dtype='float32')[c]*supported[c,None]
        lf=np.eye(len(supported_lines),dtype='float32')[line]*supported_lines[line,None]
        return self.head(features,raw[:,:,0],raw[:,:,1],torch.from_numpy(c.astype('int64')).to(device),
            torch.from_numpy(cf).to(device=device,dtype=dtype),torch.from_numpy(lf).to(device=device,dtype=dtype))

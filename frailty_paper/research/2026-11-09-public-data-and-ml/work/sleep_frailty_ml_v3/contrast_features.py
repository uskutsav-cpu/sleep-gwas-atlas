"""Verified, streamed context/contrast inputs and training-only transforms."""
from common import *

class ContrastFeatures:
    def __init__(self):
        self.data,self.original,self.manifest,self.train,self.validation,self.test=load();self.study=self.data['study'];self.nsal=int(np.sum(self.study==0));self.active=np.array([int(r['active_screen']) for r in self.manifest],dtype='int32')
        receipt=json.loads((OUT/'motif_feature_receipt.json').read_text());self.M=receipt['motif_count'];self.motifs=[]
        for stem in ['salomon','siraj']:
            r=next(r for r in receipt['datasets'] if r['study']==stem);assert r['input_sha256']==sha(WORK/(stem+'_sequences_local.npz')) and r['features_sha256']==sha(r['path']);self.motifs.append(np.load(r['path'],mmap_mode='r'))
        cr=json.loads((OUT/'reference_chromatin_feature_receipt.json').read_text());path=WORK/'reference_chromatin_features_local.npz';assert cr['feature_sha256']==sha(path)
        with np.load(path) as z:self.chromatin=z['features']
        self.dc=2*self.M+128+self.chromatin.shape[1];self.dv=self.M+384
    def raw(self,indices):
        ix=np.asarray(indices);m=np.empty((len(ix),3*self.M),dtype='float32');a=ix<self.nsal;m[a]=self.motifs[0][ix[a]];m[~a]=self.motifs[1][ix[~a]-self.nsal];h=self.original['hyena_original'][ix]
        return np.c_[m[:,:self.M],m[:,2*self.M:],h[:,:128],self.chromatin[ix]].astype('float32'),np.c_[m[:,self.M:2*self.M],h[:,128:]].astype('float32')
    def fit(self,training):
        sums=np.zeros(self.dc,dtype='float64');squares=sums.copy();count=np.zeros(self.dc,dtype='int64');vss=np.zeros(self.dv,dtype='float64');n=0
        for a in range(0,len(training),1024):
            c,v=self.raw(training[a:a+1024]);ok=np.isfinite(c);safe=np.where(ok,c,0).astype('float64');sums+=safe.sum(0);squares+=(safe*safe).sum(0);count+=ok.sum(0);assert np.isfinite(v).all() and np.min(v)>-1e-5;vss+=(v.astype('float64')**2).sum(0);n+=len(v)
        self.mean=np.divide(sums,count,out=np.zeros(self.dc),where=count>0);self.sd=np.sqrt(np.maximum(np.divide(squares,count,out=np.zeros(self.dc),where=count>0)-self.mean**2,0));self.sd[self.sd<1e-5]=1.;self.vscale=np.sqrt(vss/n);self.vscale[self.vscale<1e-5]=1.;self.mean=self.mean.astype('float32');self.sd=self.sd.astype('float32');self.vscale=self.vscale.astype('float32');return self
    def transform(self,indices):
        c,v=self.raw(indices);return np.clip((np.where(np.isfinite(c),c,self.mean)-self.mean)/self.sd,-20,20).astype('float32'),np.clip(v/self.vscale,-20,20).astype('float32')
    def load_transform(self,z):
        self.mean=z['context_mean'];self.sd=z['context_sd'];self.vscale=z['variant_scale'];return self

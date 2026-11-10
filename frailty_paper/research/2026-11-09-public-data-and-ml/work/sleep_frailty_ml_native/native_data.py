"""Pack actual single-substitution constructs without a large disk cache."""
from native_common import *
from measurement_inputs import load_data
import torch


class PackedDNA:
    def __init__(self,d,bank):
        used=np.unique(d['bank_index']);needed=set(map(int,used));n=len(bank);chunks=[None]*n
        self.length=np.zeros(n,dtype='int16');self.position=np.zeros(n,dtype='int16');self.alt=np.zeros(n,dtype='uint8')
        lut=np.full(256,255,dtype='uint8');lut[[ord(c) for c in 'AGCT']]=np.arange(4,dtype='uint8');sources=[]
        for path in [V7/'sequence_bank.csv.gz',DATA/'extra_sequence_bank.csv.gz',OLD/'OL27_development_sequence_bank.csv.gz']:
            sources.append({'path':str(path),'sha256':sha(path)})
            for r in csv_rows(path):
                i=int(r['bank_id'])
                if i not in needed:continue
                assert chunks[i] is None
                a,b=r['reference'],r['alternate'];p=int(r['changed_position']);assert len(a)==len(b)
                assert [j for j,(x,y) in enumerate(zip(a,b)) if x!=y]==[p]
                encoded=lut[np.frombuffer(a.encode(),dtype='uint8')];assert encoded.max()<4
                chunks[i]=encoded.tobytes();self.length[i]=len(a);self.position[i]=p;self.alt[i]=lut[ord(b[p])]
                restored=encoded.copy();restored[p]=self.alt[i]
                assert np.frombuffer(b'AGCT',dtype='uint8')[restored].tobytes().decode()==b
        assert all(chunks[i] is not None for i in used)
        self.offset=np.r_[0,np.cumsum(self.length.astype('int64'))]
        self.tokens=np.frombuffer(b''.join(a or b'' for a in chunks),dtype='uint8')
        assert len(self.tokens)==self.offset[-1]
        self.receipt={'physical_pairs':len(used),'reference_token_bytes':len(self.tokens),
          'length_counts':{str(k):int((self.length[used]==k).sum()) for k in np.unique(self.length[used])},
          'sources':sources,'all_ALT_sequences_reconstructed_exactly_from_source':True,
          'no_padding_or_context_invention':True}

    def batch(self,banks,device,dtype):
        length=self.length[banks];assert len(np.unique(length))==1
        ref=np.stack([self.tokens[self.offset[i]:self.offset[i+1]] for i in banks]);alt=ref.copy()
        alt[np.arange(len(banks)),self.position[banks]]=self.alt[banks]
        values=np.r_[ref,alt].astype('int64')
        return torch.nn.functional.one_hot(torch.from_numpy(values),4).permute(0,2,1).to(device=device,dtype=dtype)


def symmetric_bank_features(bank,x,extra,banks):
    raw=bank[banks]
    scalar=np.concatenate([np.c_[(raw[:,r,0]+raw[:,r,1])*.5,raw[:,r,3:16]] for r in range(3)],1)
    projected=extra[banks].astype('float32');projected[:,:,[1,3],:]=abs(projected[:,:,[1,3],:])
    return np.c_[scalar,x[banks,48:50],projected.reshape(len(banks),192)].astype('float32')


def normalization(bank,x,extra,banks):
    total=np.zeros(236);squares=np.zeros(236);n=0
    for start in range(0,len(banks),4096):
        z=symmetric_bank_features(bank,x,extra,banks[start:start+4096]).astype('float64')
        total+=z.sum(0);squares+=(z*z).sum(0);n+=len(z)
    mean=total/n;sd=np.sqrt(np.maximum(squares/n-mean*mean,0));sd[sd<1e-5]=1
    return {'mean':mean.astype('float32').tolist(),'sd':sd.astype('float32').tolist()}

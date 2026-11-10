"""Small readout features assembled from existing audited banks in RAM blocks."""
from representation_common import *
from measurement_inputs import embeddings


def bank_features(extra,bank,x,used):
    # Only feature rows used by the existing development observations are built.
    magnitude=np.full((len(bank),140),np.nan,dtype='float32')
    signed=np.full((len(bank),99),np.nan,dtype='float32')
    for start in range(0,len(used),4096):
        ids=used[start:start+4096];mean,delta=embeddings(extra,ids);raw=bank[ids]
        symmetric=np.concatenate([np.c_[(raw[:,r,0]+raw[:,r,1])*.5,raw[:,r,3:16]] for r in range(3)],axis=1)
        magnitude[ids]=np.c_[abs(delta),symmetric,x[ids,48:50]]
        signed[ids]=np.c_[delta,raw[:,:,1]-raw[:,:,0]]
    assert np.isfinite(magnitude[used]).all() and np.isfinite(signed[used]).all()
    return magnitude,signed


def training_normalization(matrix,banks,signed=False):
    total=np.zeros(matrix.shape[1]);squares=np.zeros_like(total);n=0
    for start in range(0,len(banks),4096):
        a=matrix[banks[start:start+4096]].astype('float64');total+=a.sum(0);squares+=(a*a).sum(0);n+=len(a)
    center=np.zeros_like(total) if signed else total/n
    sd=np.sqrt(np.maximum(squares/n-center*center,0));sd[sd<1e-5]=1
    return {'center':center.tolist(),'scale':sd.tolist(),'signed_antisymmetric':signed}


def design(matrix,banks,lines,supported,settings,signed=False):
    values=np.clip((matrix[banks].astype('float64')-np.asarray(settings['center']))/np.asarray(settings['scale']),-12,12)
    lf=np.eye(len(supported))[lines]*supported[lines,None]
    indices=[96,97,98] if signed else [96,110,124,97,111,125]
    interactions=(lf[:,:,None]*values[:,None,indices]).reshape(len(values),-1)
    return np.c_[values,interactions] if signed else np.c_[values,interactions,lf]

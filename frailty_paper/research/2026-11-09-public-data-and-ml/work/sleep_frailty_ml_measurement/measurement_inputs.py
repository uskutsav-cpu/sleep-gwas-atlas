"""Read prior features in blocks, retaining original data/partitions unchanged."""
from measurement_common import *
from coverage_inputs import load_inputs

def load_data():
    d,pairs,bank,x=load_inputs();receipt=json.loads((CHANNEL/'channel_feature_receipt.json').read_text());assert receipt['passed'] and sha(CHANNEL/'projected_reader_channels.npy')==receipt['features_sha256'];extra=np.load(CHANNEL/'projected_reader_channels.npy',mmap_mode='r');return d,pairs,bank,x,extra

def embeddings(extra,banks):
    values=extra[banks].astype('float32');mean=np.concatenate([values[:,:,0,:],values[:,:,2,:]],axis=2).reshape(-1,96);difference=(np.concatenate([values[:,:,1,:],values[:,:,3,:]],axis=2)/1000.).reshape(-1,96);return mean,difference

def normalization(extra,x,banks):
    total=np.zeros(96);squares=np.zeros(96);delta_squares=np.zeros(96);n=0
    for i in range(0,len(banks),4096):
        mean,delta=embeddings(extra,banks[i:i+4096]);mean=mean.astype('float64');delta=delta.astype('float64');assert np.isfinite(mean).all() and np.isfinite(delta).all();total+=mean.sum(0);squares+=(mean*mean).sum(0);delta_squares+=(delta*delta).sum(0);n+=len(mean)
    center=total/n;sd=np.sqrt(np.maximum(squares/n-center*center,0));dsd=np.sqrt(delta_squares/n);sd[sd<1e-5]=1;dsd[dsd<1e-5]=1;physical=x[banks,48:50];pc=physical.mean(0);ps=physical.std(0);ps[ps<1e-5]=1
    return {'mean_center':center.astype('float32').tolist(),'mean_sd':sd.astype('float32').tolist(),'difference_RMS':dsd.astype('float32').tolist(),'physical_center':pc.tolist(),'physical_sd':ps.tolist()}

def normalized_inputs(extra,x,banks,settings):
    mean,delta=embeddings(extra,banks);mean=np.clip((mean-np.asarray(settings['mean_center'],dtype='float32'))/np.asarray(settings['mean_sd'],dtype='float32'),-12,12);delta=np.clip(delta/np.asarray(settings['difference_RMS'],dtype='float32'),-12,12);physical=np.clip((x[banks,48:50]-np.asarray(settings['physical_center'],dtype='float32'))/np.asarray(settings['physical_sd'],dtype='float32'),-12,12)
    return mean,delta,physical

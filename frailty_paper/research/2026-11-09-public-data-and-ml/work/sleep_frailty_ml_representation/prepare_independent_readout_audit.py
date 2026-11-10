"""Export a small actual-data numerical check for independent base R."""
from representation_common import *
from features import bank_features,design
from measurement_inputs import load_data
from training_utilities import transformations


def main():
    assert_versions_all();record=json.loads((DEST/'none_diagnostic.json').read_text())
    saved=np.load(DEST/'none_readouts.npz');assert sha(DEST/'none_readouts.npz')==record['readouts_sha256']
    d,pairs,bank,x,extra=load_data();load_signed_sidecar(d);tr=d['split']==0
    known=np.isfinite(d['signed_effect'])&tr;ids=[]
    for line in np.unique(d['line'][known]):ids.extend(np.flatnonzero(known&(d['line']==line))[:8])
    ids=np.array(ids);assert len(ids)>20
    mag,sign=bank_features(extra,bank,x,np.unique(d['bank_index'][ids]))
    z=design(mag,d['bank_index'][ids],d['line'][ids],np.array(record['supported_lines']),record['magnitude_normalization'])
    sz=design(sign,d['bank_index'][ids],d['line'][ids],np.array(record['signed_supported_lines']),record['signed_normalization'],True)
    med,support,weights,centers,quality,context=transformations(d,tr,True)
    line_support=np.array(record['supported_lines']);lm=np.array(record['line_medians'])
    base=np.where(line_support[d['line'][ids]],lm[d['line'][ids]],med[d['cell'][ids]])
    beta=saved['magnitude_coefficient'];signed_beta=saved['signed_coefficient'];nl=len(lm);penalty=np.r_[np.full(len(beta)-nl,.1),np.zeros(nl)]
    residual=base+z@beta-d['effect'][ids];root=np.sqrt(residual**2+.01**2);w=weights[ids].astype('float64')
    objective=float(w@root/w.sum()+.5*penalty@(beta*beta));gradient=z.T@(w*residual/root)/w.sum()+penalty*beta
    # Saved sufficient statistics independently certify the signed ridge solve.
    station=(saved['signed_gram']+.1*np.eye(len(signed_beta)))@signed_beta-saved['signed_rhs'];assert abs(station).max()<1e-10
    for name,values in [('audit_magnitude_design.csv',z),('audit_signed_design.csv',sz)]:
        np.savetxt(DEST/name,values,delimiter=',',fmt='%.17g')
    np.savetxt(DEST/'audit_magnitude_parameters.csv',np.c_[beta,penalty,gradient],delimiter=',',fmt='%.17g')
    np.savetxt(DEST/'audit_actual_labels.csv',np.c_[base,d['effect'][ids],w,d['signed_effect'][ids]],delimiter=',',fmt='%.17g')
    np.savetxt(DEST/'audit_signed_ridge_system.csv',np.c_[saved['signed_gram']+.1*np.eye(len(signed_beta)),saved['signed_rhs'],signed_beta],delimiter=',',fmt='%.17g')
    write_json(DEST/'readout_audit_expectations.json',{'magnitude_batch_objective':objective,'magnitude_batch_rows':len(ids),
       'maximum_signed_ridge_stationarity_error':float(abs(station).max()),'gradient_is_of_smoothed_absolute_error':True,
       'actual_training_observation_indices':ids.tolist(),'reserved_outcomes_read':False,'saved_weights_changed':False})
    print('ACTUAL-DATA INDEPENDENT READOUT AUDIT INPUTS PREPARED',len(ids),flush=True)


if __name__=='__main__':main()

"""Independent fitted NumPy/gradient audit, including the appended AVI inputs."""
from final_validation_common import *
from feature_banks import avi_bank


def numpy_forward(s,g,l,c,cf,lf,avi):
    def bn(x,p):return (x-s[p+'.running_mean'])/np.sqrt(s[p+'.running_var']+1e-5)*s[p+'.weight']+s[p+'.bias']
    def silu(x):return x/(1+np.exp(np.clip(-x,-600,600)))
    def embed(x):
        x=bn(x,'reader.mapper.block.0')@s['reader.mapper.block.1.weight'][:,:,0].T+s['reader.mapper.block.1.bias']
        return silu(bn(x@s['reader.head.0.weight'].T+s['reader.head.0.bias'],'reader.head.1'))
    a,b=embed(g[:,0]),embed(g[:,1]);la,lb=embed(l[:,0]),embed(l[:,1]);w=s['cell_activity.weight'][c]+lf@s['line_activity.weight'];bias=s['cell_activity.bias'][c]+lf@s['line_activity.bias']
    ar=(a*w).sum(1)+bias;aa=(b*w).sum(1)+bias;mean=(a+b)*.5;centered=mean-mean.mean(1,keepdims=True);normalized=centered/np.sqrt((centered**2).mean(1,keepdims=True)+1e-5)
    x=np.c_[abs(b-a)*10.,normalized,abs(aa-ar)*10.,(ar+aa)*.5,np.log1p(np.sqrt((mean*mean).mean(1)+1e-12)),abs(lb-la)*10.,cf,lf]
    if avi is not None:x=np.c_[x,avi]
    h=silu(x@s['mag_hidden.weight'].T+s['mag_hidden.bias']);raw=s['mag_cell_bias'][c]+lf@s['mag_line_offset']+h@s['mag_out.weight'][0]+s['mag_out.bias'][0]
    return {'magnitude':np.logaddexp(0,raw),'signed_effect':aa-ar,'activity_reference':ar,'activity_alternate':aa}


def main():
    assert_previous_sources();with_data=np.load(DATA/'training_arrays.npz');d={k:with_data[k] for k in with_data.files};with_data.close();checks={}
    for arm in ['sequence','AVI']:
        r,models,ensemble=load_ensemble('none',arm);ids=ensemble['observation_index'][np.linspace(0,len(ensemble['observation_index'])-1,24,dtype=int)]
        g=global_bank()[d['bank_index'][ids]].astype('float64');l=local_bank()[d['bank_index'][ids]].astype('float64');c=d['cell'][ids].astype('int64');line=d['line'][ids].astype('int64')
        cf=np.eye(len(d['cell_names']))[c]*np.array(r['supported_cells'])[c,None];lf=np.eye(len(d['cell_line_names']))[line]*np.array(r['supported_lines'])[line,None]
        av=None
        if arm=='AVI':
            t=r['AVI_transform'];av=avi_bank()[d['bank_index'][ids]].astype('float64');av[:,:2]=(np.where(np.isfinite(av[:,:2]),av[:,:2],t['median'])-t['mean'])/t['sd'];av=np.clip(av,-12,12)
        model=models[0].double();s={k:v.detach().numpy().copy() for k,v in model.state_dict().items()};expected=numpy_forward(s,g,l,c,cf,lf,av)
        tensors=[torch.from_numpy(v) for v in [g[:,0],g[:,1],l[:,0],l[:,1],c,line,cf,lf]]
        def predict():return model(*tensors,torch.from_numpy(av) if av is not None else None)
        got=predict();errors={k:float(abs(got[k].detach().numpy()-expected[k]).max()) for k in expected};assert max(errors.values())<1e-7,errors
        with np.load(FINAL/'fits'/f'fixed_none_{arm}_42_predictions.npz') as z:
            stored=z['magnitude'][np.linspace(0,len(z['magnitude'])-1,24,dtype=int)]
        replay_error=float(abs(stored-expected['magnitude']).max());assert replay_error<2e-5,replay_error
        check={'NumPy_forward_maximum_errors':errors,'stored_prediction_replay_maximum_error':replay_error}
        if arm=='AVI':
            y=torch.from_numpy(d['effect'][ids].astype('float64'));loss=(got['magnitude']-y).abs().mean();loss.backward();parameter=model.mag_hidden.weight;gradient=parameter.grad[:,-3:];index=int(gradient.abs().argmax());row,col=np.unravel_index(index,gradient.shape);col+=parameter.shape[1]-3;analytic=float(parameter.grad[row,col]);assert abs(analytic)>1e-12
            original=float(parameter.detach()[row,col]);epsilon=1e-6
            with torch.no_grad():
                parameter[row,col]=original+epsilon;plus=float((predict()['magnitude']-y).abs().mean());parameter[row,col]=original-epsilon;minus=float((predict()['magnitude']-y).abs().mean());parameter[row,col]=original
            finite=(plus-minus)/(2*epsilon);relative=abs(finite-analytic)/max(abs(finite),abs(analytic),1e-12);assert relative<1e-4
            check.update(fitted_AVI_weight_nonzero=bool(np.any(s['mag_hidden.weight'][:,-3:]!=0)),AVI_gradient_relative_error=relative)
        checks[arm]=check
    write_json(FINAL/'final_fitted_network_independent_audit.json',{'passed':True,'checks':checks,'scope':'Fitted numerical/gradient correctness, not evidence of successful molecular or disease prediction.','fresh_effects_read':False})
    print('FINAL MATCHED NETWORKS INDEPENDENTLY VERIFIED',flush=True)


if __name__=='__main__':main()

"""Separate-locus probability calibration, prediction sets and abstention."""
from data import *
from predict_checkpoints import predict_registry
from scipy.optimize import minimize
def logit(p):return np.log(np.clip(p,1e-7,1-1e-7)/(1-np.clip(p,1e-7,1-1e-7)))
def fit_logistic(p,y):
    x=logit(p);y=np.asarray(y,dtype='float64');rate=np.clip(y.mean(),1e-6,1-1e-6)
    def fun(v):
        z=v[0]*x+v[1];q=1/(1+np.exp(-np.clip(z,-40,40)));loss=np.sum(np.logaddexp(0,z)-y*z)+.5*v[0]**2;grad=np.array([np.sum((q-y)*x)+v[0],np.sum(q-y)])
        return loss,grad
    result=minimize(fun,[1.,np.log(rate/(1-rate))],jac=True,method='L-BFGS-B',bounds=[(0,None),(None,None)],options={'maxiter':1000,'ftol':1e-12})
    assert result.success and np.isfinite(result.x).all()
    return {'slope':float(result.x[0]),'intercept':float(result.x[1]),'rows':len(y),'calls':int(y.sum()),'monotonic':True}
def apply_logistic(p,c):return 1/(1+np.exp(-np.clip(c['slope']*logit(p)+c['intercept'],-40,40)))
def grouped_quantile(scores,groups,alpha=.1):
    unique,inverse=np.unique(groups,return_inverse=True);maximum=np.full(len(unique),-np.inf);np.maximum.at(maximum,inverse,scores);k=int(np.ceil((len(maximum)+1)*(1-alpha)));q=1. if k>len(maximum) else float(np.sort(maximum)[k-1])
    return {'q':q,'calibration_loci':len(unique),'observations':len(scores),'quantile_order':k,'alpha':alpha,'unit':'Maximum nonconformity per locus, class and applicable cell context'}
def apply_policy(raw,cells,distance,policy):
    p=np.empty(len(raw));include0=np.empty(len(raw),dtype=bool);include1=np.empty(len(raw),dtype=bool);reason=[]
    for cell in np.unique(cells):
        use=cells==cell;cal=policy['cell_probability'].get(str(cell),policy['pooled_probability']);p[use]=apply_logistic(raw[use],cal);threshold=policy['cell_prediction_sets'].get(str(cell),policy['pooled_prediction_sets']);include0[use]=p[use]<=threshold['0']['q'];include1[use]=(1-p[use])<=threshold['1']['q']
    known=set(policy['supported_cell_groups']);fitted=set(policy.get('fitted_cell_groups',policy['supported_cell_groups']))
    for i,cell in enumerate(cells):
        if str(cell) not in fitted:reason.append('unranked_unsupported_cell');include0[i]=True;include1[i]=True
        elif str(cell) not in known:reason.append('unranked_unvalidated_cell');include0[i]=True;include1[i]=True
        elif distance[i]>policy['domain_distance'].get(str(cell),policy['pooled_domain_distance']):reason.append('unranked_sequence_domain_shift')
        elif include0[i] and include1[i]:reason.append('unranked_weak_evidence')
        elif not include0[i] and not include1[i]:reason.append('unranked_inconsistent_evidence')
        else:reason.append('ranked')
    reason=np.array(reason);ranked=reason=='ranked'
    return {'calibrated_call_probability':p,'set_includes0':include0,'set_includes1':include1,'status':reason,'ranked_score':np.where(ranked,p,np.nan)}
def main():
    assert not (OUT/'final_test_evaluation_started.json').exists();registry=json.loads((OUT/'final_fit_registry.json').read_text());assert registry['calibration_rows_fitted'] is False
    d=load();ids=np.flatnonzero(d['split']==2);cells=np.array([d['cell_names'][i] for i in d['cell'][ids]]);feature_receipt=json.loads((OUT/'fixed_features_receipt.json').read_text());bank_x=np.load(feature_receipt['feature_path'],mmap_mode='r');predpath=OUT/'calibration_raw_predictions.npz';signature=sha(OUT/'final_fit_registry.json')+'_'+feature_receipt['feature_sha256']
    if predpath.exists():
        with np.load(predpath) as z:assert str(z['signature'])==signature and np.array_equal(z['row_index'],ids);predictions={k.removeprefix('family::'):z[k] for k in z.files if k.startswith('family::')};distance=z['distance']
    else:
        predictions,members,distance=predict_registry(registry,bank_x,d['bank_index'][ids],cells);np.savez_compressed(predpath,signature=np.array(signature),row_index=ids,distance=distance,**{'family::'+k:v for k,v in predictions.items()})
    supported=set()
    for r in registry['records']:
        if r['name'] not in registry['selected_family_weights']:continue
        path=Path(r['checkpoint'])
        if path.suffix=='.npz':
            with np.load(path) as z:fit=z['fit_row_indices']
        else:
            from predict_checkpoints import load_member
            fit=load_member(r)['fit_row_indices']
        supported.update(d['cell_names'][i] for i in d['cell'][fit])
    role=np.array([int(hashlib.sha256(('V5_probability_vs_sets|'+str(g)).encode()).hexdigest()[:8],16)%2 for g in d['groups'][ids]],dtype='int8');known=np.array([c in supported for c in cells]);prob=(role==0)&known;sets=(role==1)&known;assert not set(d['groups'][ids[prob]])&set(d['groups'][ids[sets]])
    raw=predictions['primary'];truth=d['y'][ids];pooled=fit_logistic(raw[prob],truth[prob]);cell_prob={};distance_limits={}
    for cell in sorted(supported):
        use=prob&(cells==cell)
        if use.sum()>=500 and truth[use].sum()>=50 and (use.sum()-truth[use].sum())>=50:cell_prob[str(cell)]=fit_logistic(raw[use],truth[use])
        if use.sum()>=50:distance_limits[str(cell)]=float(np.quantile(distance[use],.99))
    calibrated=np.empty(len(raw))
    for cell in np.unique(cells):use=cells==cell;calibrated[use]=apply_logistic(raw[use],cell_prob.get(str(cell),pooled))
    pooled_sets={}
    for label in [0,1]:
        use=sets&(truth==label);score=calibrated[use] if label==0 else 1-calibrated[use];pooled_sets[str(label)]=grouped_quantile(score,d['groups'][ids[use]])
    cell_sets={}
    for cell in sorted(supported):
        results={}
        for label in [0,1]:
            use=sets&(cells==cell)&(truth==label)
            if len(np.unique(d['groups'][ids[use]]))<30:break
            score=calibrated[use] if label==0 else 1-calibrated[use];results[str(label)]=grouped_quantile(score,d['groups'][ids[use]])
        if len(results)==2:cell_sets[str(cell)]=results
    trusted=[];cell_calibration_inventory={}
    for cell in sorted(supported):
        a=prob&(cells==cell);b=sets&(cells==cell);aa=len(np.unique(d['groups'][ids[a]]));bb=len(np.unique(d['groups'][ids[b]]));enough=aa>=10 and bb>=10 and a.sum()>=50 and b.sum()>=50 and min(truth[a].sum(),a.sum()-truth[a].sum(),truth[b].sum(),b.sum()-truth[b].sum())>=10
        if enough:trusted.append(cell)
        cell_calibration_inventory[cell]={'probability_role_rows':int(a.sum()),'prediction_set_role_rows':int(b.sum()),'probability_role_loci':aa,'prediction_set_role_loci':bb,'eligible_for_ranked_output':bool(enough)}
    policy={'primary_registry_sha256':sha(OUT/'final_fit_registry.json'),'calibration_prediction_sha256':sha(predpath),'fitted_cell_groups':sorted(supported),'supported_cell_groups':trusted,'cell_calibration_inventory':cell_calibration_inventory,'support_rule':'Rank only fitted cell groups with >=10distinct calibration loci,>=50observations and>=10ofeachlabel in EACH separate calibration role. Other trained cells are unvalidated/unranked.','pooled_probability':pooled,'cell_probability':cell_prob,'pooled_prediction_sets':pooled_sets,'cell_prediction_sets':cell_sets,'domain_distance':distance_limits,'pooled_domain_distance':float(np.quantile(distance[prob],.99)),'probability_role_rows':int(prob.sum()),'set_role_rows':int(sets.sum()),'unsupported_calibration_rows':int((~known).sum()),'role_assignment':'SHA256 V5_probability_vs_sets|locus_group mod2,0probability,1sets','same_locus_used_for_both_roles':False,'nominal_label_coverage':.9,'limitations':'Locus maxima address correlated rows within loci, but exchangeable loci/cells and stable study distributions are not guaranteed. Nominal coverage is not proof of actual new-study coverage.','final_study_outcomes_used':False,'probability_calibration_scope':'Probability of a source-style statistical laboratory detection, not biological causality or patient disease risk.'}
    write_json(OUT/'uncertainty_policy.json',policy);result=apply_policy(raw,cells,distance,policy);np.savez_compressed(OUT/'calibration_uncertainty_predictions.npz',row_index=ids,role=role,**result)
    print('UNCERTAINTY POLICY FROZEN',policy['probability_role_rows'],policy['set_role_rows'],'supported cells',policy['supported_cell_groups'],flush=True)
if __name__=='__main__':main()

"""Source-protected error learning and worst-source scaling, before reserved outcomes."""
import sys,time,argparse,joblib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
from coverage_inputs import load_inputs
from coverage_validation_common import *
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_absolute_error

def fit_error_package(features,y,p,baseline,cell,source,groups,fit,scale,policy):
    assert fit.any() and scale.any();assert not set(groups[fit])&set(groups[scale]);critics={};strata=source[fit].astype('int64')*100+cell[fit];_,inverse,counts=np.unique(strata,return_inverse=True,return_counts=True);weight=1./counts[inverse];weight/=weight.mean();error=abs(y-p)
    for label,target in [('mean_error',error),('baseline_error',abs(y-baseline)),('quantile_error',error)]:
        params=policy['critic_parameters'].copy()
        if label=='quantile_error':params.update(loss='quantile',quantile=.9)
        critic=HistGradientBoostingRegressor(**params);critic.fit(features[fit],target[fit],sample_weight=weight);critics[label]=critic
    raw=np.maximum(critics['quantile_error'].predict(features[scale]),1e-4);ratio=error[scale]/raw;cg=groups[scale];cs=source[scale];cc=cell[scale];ce=error[scale]
    def conservative(mask):
        assert mask.any();factor=weighted_quantile(ratio[mask],cg[mask]);floor=weighted_quantile(ce[mask],cg[mask]);details=[]
        for sid in np.unique(cs[mask]):
            ix=mask&(cs==sid);ng=len(np.unique(cg[ix]))
            if ix.sum()>=policy['minimum_scale_rows'] and ng>=policy['minimum_scale_regions']:
                f=weighted_quantile(ratio[ix],cg[ix]);e=weighted_quantile(ce[ix],cg[ix]);factor=max(factor,f);floor=max(floor,e);details.append({'source':int(sid),'rows':int(ix.sum()),'regions':ng,'factor':f,'absolute_error_floor':e})
        return {'factor':max(factor,1e-4),'absolute_error_floor':max(floor,1e-4),'source_details':details,'rows':int(mask.sum()),'regions':len(np.unique(cg[mask]))}
    global_scale=conservative(np.ones(len(raw),dtype=bool));by_cell={}
    for c in np.unique(cell):
        mask=cc==c
        if mask.sum()>=policy['minimum_scale_rows'] and len(np.unique(cg[mask]))>=policy['minimum_scale_regions']:by_cell[int(c)]=conservative(mask)
        else:by_cell[int(c)]=dict(global_scale,scope='global_fallback')
    center=features[fit].mean(0);sd=features[fit].std(0);sd[sd<1e-5]=1
    return {'critics':critics,'by_cell':by_cell,'global_scale':global_scale,'feature_center':center,'feature_sd':sd,
            'fit_rows':int(fit.sum()),'scale_rows':int(scale.sum()),'fit_groups':len(np.unique(groups[fit])),'scale_groups':len(np.unique(groups[scale]))}

def uncertainty(package,features,cell,baseline,cs,ls,counts,policy):
    critics=package['critics'];raw=np.maximum(critics['quantile_error'].predict(features),1e-4);scales=[package['by_cell'].get(int(c),package['global_scale']) for c in cell]
    width=np.maximum(raw*np.array([s['factor'] for s in scales]),np.array([s['absolute_error_floor'] for s in scales]));estimated=np.maximum(critics['mean_error'].predict(features),0);base=np.maximum(critics['baseline_error'].predict(features),1e-4);z=np.max(abs((features-package['feature_center'])/package['feature_sd']),axis=1)
    accepted=cs&ls&(counts>=policy['minimum_training_studies_same_cell'])&(z<=policy['maximum_standardized_input'])&((base-estimated)/base>=policy['minimum_estimated_relative_gain'])&((base-estimated)>=policy['minimum_estimated_absolute_gain'])&(width<=np.maximum(policy['maximum_halfwidth_absolute_floor'],baseline))
    return {'raw_error_scale':raw.astype('float32'),'interval_halfwidth':width.astype('float32'),'estimated_error':estimated.astype('float32'),'estimated_baseline_error':base.astype('float32'),'maximum_standardized_error_feature':z.astype('float32'),'pretest_policy_accepted':accepted}

def summary(y,p,b,width,accepted):
    result={'rows':len(y),'MAE':float(abs(y-p).mean()),'baseline_MAE':float(abs(y-b).mean()),'coverage':float(np.mean(abs(y-p)<=width)),
            'median_halfwidth':float(np.median(width)),'accepted_rows':int(accepted.sum()),'accepted_fraction':float(accepted.mean()),'interpretation':'Development diagnostic, not untouched validation.'}
    if accepted.any():result['accepted']={'MAE':float(abs(y[accepted]-p[accepted]).mean()),'baseline_MAE':float(abs(y[accepted]-b[accepted]).mean()),'coverage':float(np.mean(abs(y[accepted]-p[accepted])<=width[accepted]))}
    return result

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--arm',choices=['sequence','AVI'],required=True);args=parser.parse_args();start=time.monotonic();assert_previous_sources()
    assert not (FINAL/'before_reserved_outcomes_lock.json').exists();assert (FINAL/(args.arm+'_fixed_execution_receipt.json')).exists()
    policy={'critic_parameters':{'max_iter':80,'max_leaf_nodes':7,'min_samples_leaf':100,'l2_regularization':10.,'early_stopping':False,'random_state':42},
            'minimum_scale_rows':100,'minimum_scale_regions':20,'minimum_training_studies_same_cell':1,'maximum_standardized_input':10.,
            'minimum_estimated_relative_gain':.05,'minimum_estimated_absolute_gain':.005,'maximum_halfwidth_absolute_floor':.2,
            'scale_method':'For each supported cell, take max of pooled and each adequately represented source 90th region-balanced residual/scale quantile; also max absolute residual 90th quantile as width floor. Sparse cells use global fallback.',
            'acceptance':'Known exact cell line and cell, within standardized input bounds, estimated error gain >=5% and >=0.005 log2 units, interval halfwidth <=max(0.2 log2 units, training-only typical magnitude). These estimates never authorize release without fresh accuracy/coverage tests.',
            'nested_split':'Other-source calibration groups from outer-source-omitted base model: SHA256 error-partition-20261108|group even error fit, odd interval scale. Outer source absent from base fitting, error fitting, and interval scaling.',
            'final_split':'Whole-source-omission predictions; original noncalibration regions fit critics, original calibration regions scale intervals. No reserved outcomes.',
            'coverage_claim':'Nominal 90% target. No universal source transfer guarantee; report actual new-study coverage and interval size.',
            'required_fresh_observed_coverage':.9,'all_withheld_or_excessively_broad_intervals_are_not_useful_success':True}
    policy_path=FINAL/'confidence_policy_lock.json'
    if policy_path.exists():assert json.loads(policy_path.read_text())==policy
    else:write_json(policy_path,policy)
    d,pairs,bank,x=load_inputs();settings=json.loads((FINAL/'final_settings_lock.json').read_text());nested={};outer_records=[]
    for source_name in settings['held_studies'][1:]:
        meta,models,p=load_ensemble(source_name,args.arm);del models;ix=p['observation_index'];features,b,cm,cs,ls,counts=confidence_inputs(p,x[d['bank_index'][ix]],d['cell'][ix],d['line'][ix],meta)
        ometa,omodels,other=load_ensemble(source_name,args.arm,True);del omodels;oi=other['observation_index'];assert np.all(d['split'][oi]==2);sid=list(d['study_names']).index(source_name);assert not np.any(d['study'][oi]==sid)
        of,ob,_,ocs,_,_=confidence_inputs(other,x[d['bank_index'][oi]],d['cell'][oi],d['line'][oi],ometa);eligible=ocs;oi=oi[eligible];of=of[eligible];ob=ob[eligible];op=other['magnitude'][eligible]
        fit=np.array([int(hashlib.sha256(('error-partition-20261108|'+str(g)).encode()).hexdigest()[:16],16)%2==0 for g in d['group'][oi]]);scale=~fit
        package=fit_error_package(of,d['effect'][oi],op,ob,d['cell'][oi],d['study'][oi],d['group'][oi],fit,scale,policy)
        u=uncertainty(package,features,d['cell'][ix],b,cs,ls,counts,policy);nested[source_name]=summary(d['effect'][ix],p['magnitude'],b,u['interval_halfwidth'],u['pretest_policy_accepted'])
        nested[source_name].update(error_fit_rows=package['fit_rows'],error_scale_rows=package['scale_rows'],outer_source_absent_from_all_three_fitting_stages=True)
        outer_records.append((meta,p,ix,features,b,cs));print('NESTED SOURCE CONFIDENCE CHECKED',args.arm,source_name,nested[source_name]['coverage'],flush=True)
    ii=np.concatenate([r[2] for r in outer_records]);ff=np.concatenate([r[3] for r in outer_records]);bb=np.concatenate([r[4] for r in outer_records]);pp=np.concatenate([r[1]['magnitude'] for r in outer_records]);supported=np.concatenate([r[5] for r in outer_records]);assert len(np.unique(ii))==len(ii)
    ii=ii[supported];ff=ff[supported];bb=bb[supported];pp=pp[supported];fit=d['split'][ii]!=2;scale=~fit
    package=fit_error_package(ff,d['effect'][ii],pp,bb,d['cell'][ii],d['study'][ii],d['group'][ii],fit,scale,policy);package.update(policy=policy,settings_lock_sha256=sha(FINAL/'final_settings_lock.json'),policy_lock_sha256=sha(policy_path),fit_observation_index_sha256=hashlib.sha256(ii[fit].tobytes()).hexdigest(),scale_observation_index_sha256=hashlib.sha256(ii[scale].tobytes()).hexdigest())
    # Direct pair residual scaling is estimated from the same original calibration regions.
    inverse=np.full(len(d['effect']),-1,dtype='int32');inverse[ii]=np.arange(len(ii));raw=np.maximum(package['critics']['quantile_error'].predict(ff),1e-4);pair_scale={}
    for name,p in pairs.items():
        p=p[np.all(inverse[p]>=0,axis=1)&np.all(d['split'][p]==2,axis=1)];a,b=p[:,0],p[:,1];ia,ib=inverse[a],inverse[b];residual=abs((d['effect'][a]-d['effect'][b])-(pp[ia]-pp[ib]));rs=raw[ia]+raw[ib];groups=d['group'][a];factor=weighted_quantile(residual/rs,groups);floor=weighted_quantile(residual,groups);details=[]
        for sid in np.unique(d['study'][a]):
            mask=d['study'][a]==sid
            if mask.sum()>=100 and len(np.unique(groups[mask]))>=20:
                f=weighted_quantile((residual/rs)[mask],groups[mask]);e=weighted_quantile(residual[mask],groups[mask]);factor=max(factor,f);floor=max(floor,e);details.append({'source':int(sid),'factor':f,'absolute_error_floor':e})
        pair_scale[name]={'factor':factor,'absolute_error_floor':floor,'pairs':len(p),'regions':len(np.unique(groups)),'source_details':details,'scope':'Direct paired-difference calibration; individual interval overlap is not used as a pair guarantee.'}
    package['paired_difference_scale']=pair_scale;path=FINAL/(args.arm+'_confidence.joblib');joblib.dump(package,path,compress=3)
    write_json(FINAL/(args.arm+'_confidence_assessment.json'),{'nested_outer_source_development':nested,'final_fit_rows':package['fit_rows'],'final_scale_rows':package['scale_rows'],
              'global_scale':package['global_scale'],'by_cell':package['by_cell'],'paired_difference_scale':pair_scale,'policy_sha256':sha(policy_path),'confidence_package_sha256':sha(path),
              'new_test_outcomes_used':False,'useful_cross_experiment_prediction_established':False,'elapsed_seconds':time.monotonic()-start})
    print('CONFIDENCE PACKAGE COMPLETE',args.arm,flush=True)

if __name__=='__main__':main()

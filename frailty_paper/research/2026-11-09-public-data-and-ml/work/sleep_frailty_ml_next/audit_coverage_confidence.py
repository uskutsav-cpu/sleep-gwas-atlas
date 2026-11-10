"""Reconstruct outer exclusions and final interval scales with separate quantile code."""
import sys,argparse,time,joblib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
from coverage_inputs import load_inputs
from coverage_validation_common import DEST,FINAL,load_ensemble,confidence_inputs

def independently_weighted_quantile(values,groups):
    counts=collections.Counter(map(int,groups));entries=sorted((float(value),1./counts[int(group)]) for value,group in zip(values,groups));target=.9*sum(w for value,w in entries);running=0.
    for value,weight in entries:
        running+=weight
        if running>=target:return value
    return entries[-1][0]

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--arm',choices=['sequence','AVI'],required=True);args=parser.parse_args();start=time.monotonic();assert_previous_sources();d,pairs,bank,x=load_inputs();settings=json.loads((FINAL/'final_settings_lock.json').read_text());assessment=json.loads((FINAL/(args.arm+'_confidence_assessment.json')).read_text());package=joblib.load(FINAL/(args.arm+'_confidence.joblib'));all_ids=[];all_f=[];all_p=[];outer=[]
    assert sha(FINAL/(args.arm+'_confidence.joblib'))==assessment['confidence_package_sha256'];assert package['policy_lock_sha256']==sha(FINAL/'confidence_policy_lock.json')
    for source in settings['held_studies'][1:]:
        sid=list(d['study_names']).index(source);outer_groups=set(d['group'][d['study']==sid]);train=(d['split']!=2)&(d['study']!=sid)&~np.isin(d['group'],list(outer_groups));tr=np.flatnonzero(train)
        meta,models,p=load_ensemble(source,args.arm);del models;assert meta['training_index_sha256']==hashlib.sha256(tr.tobytes()).hexdigest();assert not np.any(d['study'][tr]==sid)
        ometa,models,op=load_ensemble(source,args.arm,True);del models;oi=op['observation_index'];assert np.all(d['split'][oi]==2) and not np.any(d['study'][oi]==sid);assert not (set(d['group'][oi])&outer_groups or set(d['group'][oi])&set(d['group'][tr]))
        _,_,_,ocs,_,_=confidence_inputs(op,x[d['bank_index'][oi]],d['cell'][oi],d['line'][oi],ometa);oi=oi[ocs]
        fit=np.array([int(hashlib.sha256(('error-partition-20261108|'+str(g)).encode()).hexdigest()[:16],16)%2==0 for g in d['group'][oi]]);scale=~fit
        assert not set(d['group'][oi[fit]])&set(d['group'][oi[scale]]);assert assessment['nested_outer_source_development'][source]['error_fit_rows']==int(fit.sum());assert assessment['nested_outer_source_development'][source]['error_scale_rows']==int(scale.sum())
        outer.append({'source':source,'base_train_rows':len(tr),'error_fit_rows':int(fit.sum()),'interval_scale_rows':int(scale.sum()),'outer_source_present_in_any_stage':False,
                      'error_fit_index_sha256':hashlib.sha256(oi[fit].tobytes()).hexdigest(),'scale_index_sha256':hashlib.sha256(oi[scale].tobytes()).hexdigest()})
        ix=p['observation_index'];features,_,_,cs,_,_=confidence_inputs(p,x[d['bank_index'][ix]],d['cell'][ix],d['line'][ix],meta);all_ids.append(ix[cs]);all_f.append(features[cs]);all_p.append(p['magnitude'][cs])
    ii=np.concatenate(all_ids);features=np.concatenate(all_f);pred=np.concatenate(all_p);fit=d['split'][ii]!=2;scale=~fit
    assert package['fit_observation_index_sha256']==hashlib.sha256(ii[fit].tobytes()).hexdigest();assert package['scale_observation_index_sha256']==hashlib.sha256(ii[scale].tobytes()).hexdigest();assert not set(d['group'][ii[fit]])&set(d['group'][ii[scale]])
    residual=abs(d['effect'][ii[scale]]-pred[scale]);raw=np.maximum(package['critics']['quantile_error'].predict(features[scale]),1e-4);ratio=residual/raw;cg=d['group'][ii[scale]];cs=d['study'][ii[scale]];cc=d['cell'][ii[scale]];max_error=0.
    def check(mask,stored):
        nonlocal max_error
        f=independently_weighted_quantile(ratio[mask],cg[mask]);e=independently_weighted_quantile(residual[mask],cg[mask])
        for sid in np.unique(cs[mask]):
            ix=mask&(cs==sid)
            if ix.sum()>=100 and len(np.unique(cg[ix]))>=20:f=max(f,independently_weighted_quantile(ratio[ix],cg[ix]));e=max(e,independently_weighted_quantile(residual[ix],cg[ix]))
        error=max(abs(f-stored['factor']),abs(e-stored['absolute_error_floor']));max_error=max(max_error,error);assert error<1e-6,error
    check(np.ones(len(residual),dtype=bool),package['global_scale'])
    for c,stored in package['by_cell'].items():
        mask=cc==c
        if stored.get('scope')=='global_fallback':assert stored['factor']==package['global_scale']['factor'] and stored['absolute_error_floor']==package['global_scale']['absolute_error_floor']
        else:check(mask,stored)
    result={'passed':True,'outer_source_exclusions':outer,'final_error_fit_scale_regions_disjoint':True,'independent_source_worst_quantile_max_error':max_error,
            'confidence_package_sha256':sha(FINAL/(args.arm+'_confidence.joblib')),'code_sha256':sha(Path(__file__)),
            'does_not_establish_useful_transfer':True,'reserved_outcomes_used':False,'elapsed_seconds':time.monotonic()-start}
    write_json(FINAL/(args.arm+'_confidence_independent_audit.json'),result);assert_previous_sources();print('CONFIDENCE INPUTS AND QUANTILES INDEPENDENTLY VERIFIED',args.arm,flush=True)

if __name__=='__main__':main()

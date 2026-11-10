"""Calibrate from revealed studies, then lock predictions before fresh outcomes."""
from final_validation_common import *
from sklearn.ensemble import HistGradientBoostingRegressor
import joblib,sqlite3,time


def direct_pair_calibration(d,pairs,indices,pred,scale):
    inverse=np.full(len(d['effect']),-1,dtype='int32');inverse[indices]=np.arange(len(indices));results={}
    for label,p in pairs.items():
        keep=np.all(inverse[p]>=0,axis=1)&np.all(d['split'][p]==2,axis=1);p=p[keep]
        assert len(p),(label,'No source-omission calibration pairs')
        a,b=p[:,0],p[:,1];ia,ib=inverse[a],inverse[b]
        residual=abs((d['effect'][a]-d['effect'][b])-(pred[ia]-pred[ib]));rawscale=np.maximum(scale[ia]+scale[ib],1e-4)
        groups=d['group'][a];assert np.array_equal(groups,d['group'][b]);ratio=residual/rawscale
        keys=np.array(['|'.join(sorted([str(d['cell_names'][d['cell'][x]]),str(d['cell_names'][d['cell'][y]])])) for x,y in p])
        global_factor=weighted_quantile(ratio,groups);by={}
        for key in np.unique(keys):
            ix=keys==key;enough=ix.sum()>=100 and len(np.unique(groups[ix]))>=20
            by[str(key)]={'factor':weighted_quantile(ratio[ix],groups[ix]) if enough else global_factor,'scope':'exact_cell_pair' if enough else 'global_pair_type','pairs':int(ix.sum()),'regions':len(np.unique(groups[ix]))}
        results[label]={'global_factor':global_factor,'by_cell_pair':by,'pairs':len(p),'regions':len(np.unique(groups)),
                        'method':'Direct paired-difference residual calibration, not a guarantee derived from individual intervals.'}
    return results


def main():
    start=time.monotonic();assert_previous_sources();assert not (FINAL/'fresh_outcome_evaluation.json').exists()
    receipt=json.loads((FINAL/'fixed_refit_execution_receipt.json').read_text());assert receipt['count']==48 and receipt['all_matched_arm_shared_weights_rows_and_epochs_pass']
    assert not (FINAL/'before_fresh_outcomes_lock.json').exists()
    with np.load(DATA/'training_arrays.npz') as z:d={k:z[k] for k in z.files}
    with np.load(DATA/'measured_pairs.npz') as z:pairs={k:z[k] for k in z.files}
    g,l=global_bank(),local_bank();used=np.unique(d['bank_index']);static=np.full((g.shape[0],12),np.nan,dtype='float32')
    for start_at in range(0,len(used),4096):
        ids=used[start_at:start_at+4096];static[ids]=static_features(g[ids],l[ids])
    np.save(FINAL/'confidence_static_features.npy',static)
    rows=list(csv_rows(V8/'blind_inference_row_metadata.csv'));supported_ids=np.array([i for i,r in enumerate(rows) if not r['unscorable_reason']],dtype='int64');blind_rows=[rows[i] for i in supported_ids]
    features=np.load(V8/'blind_inference_pretrained_features.npy');banks=np.array([int(r['bank_index']) for r in blind_rows]);blind_features=features[banks];blind_static=static_features(blind_features[:,:,0],blind_features[:,:,1])
    freshavi={r['variant_GRCh38']:r for r in csv_rows(V8/'fresh_AVI_comparison_inputs.csv')}
    db=sqlite3.connect((V7/'avi_exact_scores.sqlite').as_uri()+'?mode=ro',uri=True);blind_avi=[]
    for r in blind_rows:
        if r['role']=='fresh_test':a=freshavi[r['variant_GRCh38']];blind_avi.append([float(a['raw']) if a['raw'] else np.nan,float(a['phred']) if a['phred'] else np.nan,int(a['available'])])
        else:
            a=db.execute('SELECT raw,phred,available FROM scores WHERE variant=?',(r['variant_GRCh38'],)).fetchone();blind_avi.append([a[0],a[1],a[2]] if a is not None and a[2] else [np.nan,np.nan,0])
    db.close();blind_avi=np.array(blind_avi,dtype='float32');summaries={};lockfiles={}
    settings=json.loads((FINAL/'final_settings_lock.json').read_text())
    for arm in ['sequence','AVI']:
        xs=[];ys=[];baselines=[];indices=[];magnitudes=[];source_ids=[]
        for study in settings['held_studies'][1:]:
            metadata,models,pred=load_ensemble(study,arm);ix=pred['observation_index'];x,b,_,cs,_,_,names=confidence_features(pred,static[d['bank_index'][ix]],d['cell'][ix],d['line'][ix],metadata)
            eligible=cs;xs.append(x[eligible]);ys.append(d['effect'][ix][eligible]);baselines.append(b[eligible]);indices.append(ix[eligible]);magnitudes.append(pred['magnitude'][eligible]);source_ids.append(d['study'][ix][eligible]);del models
        x=np.concatenate(xs);y=np.concatenate(ys);b=np.concatenate(baselines);ix=np.concatenate(indices);p=np.concatenate(magnitudes);sources=np.concatenate(source_ids)
        assert len(np.unique(ix))==len(ix);train=d['split'][ix]==0;cal=d['split'][ix]==2;assert train.any() and cal.any()
        assert not set(d['group'][ix[train]])&set(d['group'][ix[cal]])
        strata=sources[train].astype('int64')*100+d['cell'][ix[train]];_,inv,counts=np.unique(strata,return_inverse=True,return_counts=True);weights=1./counts[inv];weights/=weights.mean()
        error=abs(y-p);base_error=abs(y-b);critics={}
        for label,target in [('mean_error',error),('baseline_error',base_error),('quantile_error',error)]:
            params={'max_iter':100,'max_leaf_nodes':7,'min_samples_leaf':100,'l2_regularization':10.,'early_stopping':False,'random_state':42}
            if label=='quantile_error':params.update(loss='quantile',quantile=.9)
            critic=HistGradientBoostingRegressor(**params);critic.fit(x[train],target[train],sample_weight=weights);critics[label]=critic
        rawscale=np.maximum(critics['quantile_error'].predict(x),1e-4)
        calibration=calibration_factors(error[cal]/rawscale[cal],d['cell'][ix[cal]],d['group'][ix[cal]],d['cell_names'])
        paircal=direct_pair_calibration(d,pairs,ix,p,rawscale)
        center=x[train].mean(0);sd=x[train].std(0);sd[sd<1e-5]=1
        package={'critics':critics,'feature_names':names,'feature_center':center,'feature_sd':sd,'calibration':calibration,'pair_calibration':paircal,
                 'training_observation_index_sha256':hashlib.sha256(ix[train].tobytes()).hexdigest(),'calibration_observation_index_sha256':hashlib.sha256(ix[cal].tobytes()).hexdigest(),
                 'training_weighting':'Equal total weight per revealed study/cell stratum. No fresh outcomes.',
                 'calibration_scope':'Original calibration roles of whole-source-omission predictions. These earlier revealed studies are not new independent validation.'}
        package_path=FINAL/f'{arm}_confidence.joblib';joblib.dump(package,package_path,compress=3)
        write_json(FINAL/f'{arm}_confidence_calibration.json',{'individual':calibration,'paired':paircal,'training_rows':int(train.sum()),'calibration_rows':int(cal.sum()),'feature_names':names,'training_weighting':package['training_weighting']})
        metadata,models,mainpred=load_ensemble('none',arm);ci=mainpred['observation_index'];cx,cb,cmed,cs,ls,counts,_=confidence_features(mainpred,static[d['bank_index'][ci]],d['cell'][ci],d['line'][ci],metadata)
        def uncertainty(x,cell,cs,ls,counts,conditions,activity_known):
            scale=np.maximum(critics['quantile_error'].predict(x),1e-4);factor=np.array([calibration['by_cell'][str(d['cell_names'][c])] for c in cell]);width=scale*factor
            e=np.maximum(critics['mean_error'].predict(x),0);be=np.maximum(critics['baseline_error'].predict(x),1e-4);zmax=np.max(abs((x-center)/sd),axis=1)
            accepted=cs&ls&(counts>=1)&(zmax<=10)&((be-e)/be>=.05)&((be-e)>=.005)&np.isin(conditions,['baseline','untreated','unreported'])&activity_known
            return {'raw_error_scale':scale,'interval_halfwidth':width,'estimated_error':e,'estimated_baseline_error':be,'maximum_standardized_error_feature':zmax,'pretest_policy_accepted':accepted}
        cal_unc=uncertainty(cx,d['cell'][ci],cs,ls,counts,np.repeat('baseline',len(ci)),np.ones(len(ci),dtype=bool))
        # Calibration performance here is descriptive, since its outcomes informed
        # interval scaling through source-omitted predictions.
        np.savez_compressed(FINAL/f'{arm}_main_calibration_predictions.npz',**mainpred,baseline=cb,cell_baseline=cmed,**cal_unc)
        cell=np.array([list(d['cell_names']).index(r['cell_group']) for r in blind_rows],dtype='int64');line=np.array([list(d['cell_line_names']).index(r['cell_line']) for r in blind_rows],dtype='int64')
        bp=predict_blind(models,metadata,blind_features,cell,line,blind_avi);bx,bb,bcm,bcs,bls,bcnt,_=confidence_features(bp,blind_static,cell,line,metadata)
        known=np.array([r['source_assay_activity_known']=='True' for r in blind_rows]);conditions=np.array([r['condition'] for r in blind_rows]);bu=uncertainty(bx,cell,bcs,bls,bcnt,conditions,known)
        np.savez_compressed(FINAL/f'{arm}_blinded_predictions.npz',row_index=supported_ids,**bp,baseline=bb,cell_baseline=bcm,**bu)
        summaries[arm]={'OOF_confidence_training_rows':int(train.sum()),'OOF_confidence_calibration_rows':int(cal.sum()),'blind_supported_rows':len(blind_rows),'fresh_rows':sum(r['role']=='fresh_test' for r in blind_rows),
                        'fresh_pretest_policy_accepted':int(bu['pretest_policy_accepted'][np.array([r['role']=='fresh_test' for r in blind_rows])].sum()),
                        'project_pretest_policy_accepted':int(bu['pretest_policy_accepted'][np.array([r['role']=='project' for r in blind_rows])].sum()),
                        'project_restriction':'Native genomic contexts lack verified source-assay activity; active-only fresh validation cannot authorize project answers.'}
        print('CONFIDENCE AND BLIND PREDICTIONS LOCKABLE',arm,flush=True)
    write_json(FINAL/'confidence_and_blind_prediction_receipt.json',{'arms':summaries,'fresh_effects_read':False,'elapsed_seconds':time.monotonic()-start,'environment':{'python':sys.version,'numpy':np.__version__,'torch':torch.__version__}})
    protocol=json.loads((V8/'outputs/V8_external_validation_protocol.json').read_text());protocol['version']='V9 coverage-amended, before fresh outcome comparison';protocol['V9_baseline_amendment']='Also compare against stronger training-only cell/culture-line median and require both baselines to pass.'
    protocol['V9_known_line_required']=True;protocol['V9_confidence_training_weighting']='Equal study/cell stratum totals; source-omission original training roles only.'
    protocol['V9_minimum_local_confirmation']={'pairs':30,'independent_1Mb_windows':20,'reason':'Avoid treating a handful of local groups as confirmed preservation of nearby ranking.'}
    protocol['V9_local_ranking_gate']='Require both nonzero-difference ordering accuracy and paired-difference Spearman to drop by at most 0.01, in addition to sufficient local groups.'
    protocol['V9_project_release_scope']='Not authorized by this assay-active-only experiment; diagnostic predictions remain separate from answers.'
    write_json(FINAL/'external_validation_protocol.json',protocol)
    for path in FINAL.rglob('*'):
        if path.is_file() and (path.suffix in ['.joblib','.npz','.safetensors'] or path.name in ['final_settings_lock.json','external_validation_protocol.json','confidence_and_blind_prediction_receipt.json']):lockfiles[str(path)]=sha(path)
    for path in WORK.glob('*.py'):lockfiles[str(path)]=sha(path)
    for path in [V8/'blind_inference_row_metadata.csv',V8/'blind_inference_pretrained_features.npy',V8/'fresh_AVI_comparison_inputs.csv']:lockfiles[str(path)]=sha(path)
    for path in [DATA/'training_arrays.npz',DATA/'measured_pairs.npz',V8/'sources/Matoba_journal_Supplementary_Table3.xlsx',V8/'reserved_matoba_physical_design.csv',ROOT/'work/sleep_frailty_ml_v7/extract_neuroimmune_sources.py',FINAL/'final_fitted_network_independent_audit.json']:
        lockfiles[str(path)]=sha(path)
    write_json(FINAL/'before_fresh_outcomes_lock.json',{'files':lockfiles,'fresh_outcomes_joined_or_model_performance_calculated':False,'future_changes':'Cannot claim this same study untouched for a subsequently changed model.'})
    print('V9 BEFORE-FRESH-OUTCOME LOCK COMPLETE',flush=True)


if __name__=='__main__':main()

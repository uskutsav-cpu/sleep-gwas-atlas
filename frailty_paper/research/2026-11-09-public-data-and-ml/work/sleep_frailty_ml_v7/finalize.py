"""Development-only selection, refits, study-omitted calibration and blind predictions."""
import argparse,subprocess,joblib,gzip,sys
from baselines import *

def select(d):
    path=OUT/'development_selection.json'
    if path.exists():return json.loads(path.read_text())
    assert not (OUT/'reserved_evaluation_started.json').exists();plan=freeze_plan();fit=np.flatnonzero(d['split']==0);va=np.flatnonzero(d['split']==1);families={};pred={}
    for family in list(plan['recipes'])+['ridge_100','ridge_1000','ridge_10000','tree_3','tree_5']:
        members=[];parts=[]
        for seed in plan['seeds'] if family in plan['recipes'] else [None]:
            name=f'development_{family}'+(f'_{seed}' if seed is not None else '');r=json.loads((OUT/(name+'_record.json')).read_text());assert sha(r['checkpoint'])==r['checkpoint_sha256'];members.append(r)
            with np.load(MODELS/(name+'_validation.npz')) as z:assert np.array_equal(z['row_index'],va);parts.append(z['prediction'])
        pred[family]=np.mean(parts,axis=0);families[family]={'members':members,'metrics':evaluate(d,fit,va,pred[family]),'eligible':not plan['recipes'].get(family,{}).get('no_cell',False)}
    eligible=sorted([f for f in families if families[f]['eligible']],key=lambda f:families[f]['metrics']['selection_score'],reverse=True);candidates=[{'name':f,'weights':{f:1.},'metrics':families[f]['metrics']} for f in eligible]
    for weight in [.5,.75]:
        a,b=eligible[:2];name=f'top2_mix_{int(100*weight)}';pred[name]=weight*pred[a]+(1-weight)*pred[b];candidates.append({'name':name,'weights':{a:weight,b:1-weight},'metrics':evaluate(d,fit,va,pred[name])})
    candidates.sort(key=lambda r:r['metrics']['selection_score'],reverse=True);report={'selected':candidates[0],'candidates':candidates,'families':families,'plan_sha256':sha(OUT/'model_search_plan.json'),'final_study_performance_seen':False,'scope':'Held-out development loci; source libraries have been used in development, not independent whole-study evidence.'};write_json(path,report);np.savez_compressed(OUT/'development_family_predictions.npz',row_index=va,**pred);print('DEVELOPMENT SELECTION',candidates[0]['name'],candidates[0]['metrics']['selection_score'],flush=True);return report

def fit_ensemble(d,ids,selection,stage,avi=None,omit=None):
    records=[]
    for family,weight in selection['selected']['weights'].items():
        members=selection['families'][family]['members']
        for old in members:
            if family in freeze_plan()['recipes']:r=fit_nn(d,ids,np.array([],dtype=int),family,old['seed'],stage,old['selected_epoch'],avi,omit)
            else:r=fit_classical(d,ids,np.array([],dtype=int),family,stage,avi,omit)
            records.append({'weight':weight/len(members),'record':r})
    return records

def predict_ensemble(d,ids,records,avi=None):
    result=np.zeros(len(ids),dtype='float64')
    for member in records:
        r=dict(member['record'])
        for key in ['checkpoint','transform']:
            path=Path(r[key]);r[key]=str(path if path.is_absolute() else ROOT/path)
        assert sha(r['checkpoint'])==r['checkpoint_sha256'];t=transform_load(Path(r['transform']));family=r['recipe']
        if family in freeze_plan()['recipes']:
            with np.load(r['checkpoint']) as z:p={k:z[k] for k in z.files}
            v=predict(d,ids,t,p,freeze_plan()['recipes'][family]['no_cell'],avi)
        else:v=predict_classical(d,ids,t,joblib.load(r['checkpoint']),avi)
        result+=member['weight']*v
    return result

def calibrate(d,selection,arm,avi=None):
    path=OUT/(arm+'_calibration.json')
    if path.exists():return json.loads(path.read_text())
    fit=np.flatnonzero(d['split']<2);cal=np.flatnonzero((d['split']==2)&(d['unit']==0));rows=[];folds=[]
    for study in np.unique(d['study'][cal]):
        held=cal[d['study'][cal]==study];train=fit[d['study'][fit]!=study];members=fit_ensemble(d,train,selection,'study_omitted',avi,int(study));p=predict_ensemble(d,held,members,avi);folds.append({'omitted_study':d['study_names'][int(study)],'members':members,'n':len(held)})
        for i,prediction in zip(held,p):rows.append((int(i),int(d['group'][i]),int(d['cell'][i]),int(study),float(d['effect'][i]),float(prediction)))
    a=np.array(rows);group=a[:,1].astype('int64');role_a=np.array([int(hashlib.sha256(('V7_uncertainty_'+str(g)).encode()).hexdigest()[:8],16)%2==0 for g in group]);err=abs(a[:,4]-a[:,5]);global_scale=max(float(np.median(err[role_a])),.001);scales={}
    for cell in np.unique(a[:,2]):
        use=role_a&(a[:,2]==cell);scales[int(cell)]=max(float(np.median(err[use])),.001) if use.sum()>=30 else global_scale
    scale=np.array([scales[int(c)] for c in a[:,2]]);standardized=err/scale;quantiles={}
    for study in np.unique(a[:,3]):
        ix=(~role_a)&(a[:,3]==study);n=int(ix.sum());q=float(np.quantile(standardized[ix],min(1.,np.ceil((n+1)*.9)/n),method='higher'));quantiles[d['study_names'][int(study)]]={'n':n,'q':q}
    worst=max(x['q'] for x in quantiles.values());max_width=2*float(np.quantile(d['effect'][fit[d['unit'][fit]==0]],.9));support={}
    for c in range(len(d['cell_names'])):support[c]=len(np.unique(d['study'][fit[(d['unit'][fit]==0)&(d['cell'][fit]==c)]]))
    out=OUT/(arm+'_study_omitted_residuals.csv.gz')
    with gzip.open(out,'wt') as f:
        writer=csv.writer(f);writer.writerow(['row_index','locus_group','cell_code','omitted_study_code','measured_absolute_log2_effect','predicted_absolute_log2_effect','scale_role_A']);writer.writerows([list(r)+[int(flag)] for r,flag in zip(rows,role_a)])
    report={'arm':arm,'folds':folds,'n':len(rows),'cell_error_scale':scales,'global_error_scale':global_scale,'source_quantiles_B':quantiles,'worst_study_standardized_q90':worst,'maximum_interval_width':max_width,'cell_training_study_support':support,'within_crossfit_B_coverage':float(np.mean(err[~role_a]<=scale[~role_a]*worst)),'median_interval_halfwidth':float(np.median(scale*worst)),'scale_and_quantile_roles_share_loci':bool(set(group[role_a])&set(group[~role_a])),'residual_file':str(out),'residual_sha256':sha(out),'scope':'Study-omitted calibration diagnostic, not an independent test or a cross-study90%coverage guarantee. Selection family was chosen using development studies.'};assert not report['scale_and_quantile_roles_share_loci'];write_json(path,report);print('CALIBRATION COMPLETE',arm,report['median_interval_halfwidth'],flush=True);return report

def reserved_data(d):
    rows=list(csv.DictReader((DATA/'reserved_physical_design.csv').open()));x=np.load(DATA/'reserved_feature_bank.npy',mmap_mode='r');n=len(rows);test=dict(d);test['x']=x;test['n']=2*n;test['bank_index']=np.tile(np.arange(n,dtype='int32'),2);test['unit']=np.zeros(2*n,dtype='int8');test['line']=np.full(2*n,-1,dtype='int16');test['cell']=np.repeat([d['cell_names'].index('neural_progenitor'),d['cell_names'].index('kidney_derived')],n).astype('int16');test['effect']=np.zeros(2*n);test['split']=np.full(2*n,3,dtype='int8');return test,rows

def blind_predictions(d,records,calibration,arm,avi=None):
    test,rows=reserved_data(d);ids=np.arange(test['n']);p=predict_ensemble(test,ids,records,avi);member=records[0]['record'];t=transform_load(Path(member['transform']));maxz=[]
    for start in range(0,len(rows),1024):x=np.asarray(test['x'][start:start+1024]);maxz.extend(np.max(abs((x-t['mean'])/t['sd']),axis=1).tolist())
    maxz=np.tile(maxz,2);cal=calibration;out=OUT/(arm+'_reserved_blind_predictions.csv');scales=cal['cell_error_scale'];support=cal['cell_training_study_support']
    with out.open('w') as f:
        writer=csv.DictWriter(f,fieldnames=['design_index','variant_GRCh38','cell','predicted_absolute_log2_effect','lower90','upper90','rankable','unranked_reason','max_abs_standardized_feature','AVI_available']);writer.writeheader()
        for i,value in enumerate(p):
            cell=int(test['cell'][i]);s=scales.get(str(cell),scales.get(cell,cal['global_error_scale']));nsource=support.get(str(cell),support.get(cell,0));half=s*cal['worst_study_standardized_q90'];lo=max(0.,value-half);hi=value+half;reasons=[]
            if nsource<2:reasons.append('cell_group_supported_by_fewer_than_two_training_studies')
            if maxz[i]>12.:reasons.append('sequence_features_outside_frozen_domain_limit')
            if hi-lo>cal['maximum_interval_width']:reasons.append('interval_too_wide')
            writer.writerow({'design_index':i%len(rows),'variant_GRCh38':rows[i%len(rows)]['variant_GRCh38'],'cell':['HNP','HEK293'][i//len(rows)],'predicted_absolute_log2_effect':value,'lower90':lo,'upper90':hi,'rankable':int(not reasons),'unranked_reason':';'.join(reasons),'max_abs_standardized_feature':maxz[i],'AVI_available':int(avi[i%len(rows),2]) if avi is not None else ''})
    return {'file':str(out),'sha256':sha(out),'rows':len(p),'outcomes_loaded':False}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--wait-for-pid',type=int);args=parser.parse_args()
    if args.wait_for_pid:
        while subprocess.run(['ps','-p',str(args.wait_for_pid)],capture_output=True).returncode==0:time.sleep(10)
    assert not (OUT/'reserved_evaluation_started.json').exists()
    # A separate process releases memory before the neural/ensemble steps.
    subprocess.run([sys.executable,str(WORK/'baselines.py')],check=True,env=environment());d=load();selection=select(d);fit=np.flatnonzero(d['split']<2);mx.metal.set_memory_limit(1024*1024**2);mx.metal.set_cache_limit(64*1024**2)
    base=fit_ensemble(d,fit,selection,'final');base_cal=calibrate(d,selection,'sequence_only');base_pred=blind_predictions(d,base,base_cal,'sequence_only');mx.metal.clear_cache()
    while not (OUT/'AVI_collection_receipt.json').exists():print('WAITING FOR EXACT PUBLIC AVI INPUTS',flush=True);time.sleep(10)
    avi=np.load(DATA/'avi_bank.npy',mmap_mode='r');ar=json.loads((OUT/'AVI_collection_receipt.json').read_text());assert sha(DATA/'avi_bank.npy')==ar['outputs']['avi_bank.npy'];added=fit_ensemble(d,fit,selection,'final',avi);added_cal=calibrate(d,selection,'AVI_added',avi);added_pred=blind_predictions(d,added,added_cal,'AVI_added',np.load(DATA/'avi_reserved.npy'))
    lock={'sequence_only':{'records':base,'calibration_sha256':sha(OUT/'sequence_only_calibration.json'),'blind_predictions':base_pred},'AVI_added':{'records':added,'calibration_sha256':sha(OUT/'AVI_added_calibration.json'),'blind_predictions':added_pred},'development_selection_sha256':sha(OUT/'development_selection.json'),'plan_sha256':sha(OUT/'model_search_plan.json'),'training_and_feature_code_sha256':{name:sha(WORK/name) for name in ['common.py','data.py','sequence_features.py','network.py','train.py','metrics.py','baselines.py','finalize.py','prepare_development.py','pack_data.py']},'reserved_feature_receipt_sha256':sha(OUT/'reserved_feature_receipt.json'),'reserved_effect_file_sha256':sha(DATA/'reserved_original_effects.csv.gz'),'previous_V6_lock':assert_previous_intact(),'all_predictions_locked_before_test_performance':True,'final_study_performance_seen':False};write_json(OUT/'final_model_and_prediction_lock.json',lock);print('BOTH ARMS AND PREDICTIONS LOCKED BEFORE FINAL TEST',flush=True)
if __name__=='__main__':main()

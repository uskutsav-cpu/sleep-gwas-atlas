"""Single post-lock join to original coefficients; retain unfavorable results."""
import sys,time,math,joblib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
from coverage_validation_common import DEST,FINAL
from audit_and_lock_library_measurements import write_gzip_csv
from scipy.stats import spearmanr
SCREEN=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen'

def correlation(a,b):
    return float(spearmanr(a,b).statistic) if len(a)>2 and np.ptp(a)>0 and np.ptp(b)>0 else None

def bounds(values):
    a=np.asarray([v for v in values if v is not None and np.isfinite(v)])
    return np.quantile(a,[.025,.975]).tolist() if len(a)>=100 else None

def macro_error(y,p,strata):return float(np.mean([abs(y[strata==s]-p[strata==s]).mean() for s in np.unique(strata)]))

def macro_correlation(y,p,strata):
    values=[correlation(y[strata==s],p[strata==s]) for s in np.unique(strata)];values=[v for v in values if v is not None];return float(np.mean(values)) if values else None

def primary_metrics(y,p,line_base,cell_base,groups,strata,repetitions,seed,comparison=None):
    error=macro_error(y,p,strata);lb=macro_error(y,line_base,strata);cb=macro_error(y,cell_base,strata);rho=macro_correlation(y,p,strata);clusters=np.unique(groups);parts=[np.flatnonzero(groups==g) for g in clusters];rng=np.random.default_rng(seed);gains_line=[];gains_cell=[];rhos=[];compare=[]
    for iteration in range(repetitions):
        ids=np.concatenate([parts[j] for j in rng.integers(len(parts),size=len(parts))]);e=macro_error(y[ids],p[ids],strata[ids]);l=macro_error(y[ids],line_base[ids],strata[ids]);c=macro_error(y[ids],cell_base[ids],strata[ids]);gains_line.append(1-e/l if l else None);gains_cell.append(1-e/c if c else None);rhos.append(macro_correlation(y[ids],p[ids],strata[ids]))
        if comparison is not None:
            other=macro_error(y[ids],comparison[ids],strata[ids]);compare.append(1-e/other if other else None)
    result={'rows':len(y),'genomic_1Mb_windows':len(clusters),'macro_MAE':error,'macro_training_line_or_cell_median_MAE':lb,'macro_training_cell_median_MAE':cb,
            'relative_MAE_gain_vs_line_median':1-error/lb,'relative_MAE_gain_vs_cell_median':1-error/cb,'gain_vs_line_95CI':bounds(gains_line),'gain_vs_cell_95CI':bounds(gains_cell),
            'macro_within_cell_line_Spearman':rho,'macro_within_line_Spearman_95CI':bounds(rhos),
            'cell_lines':{str(s):{'rows':int((strata==s).sum()),'MAE':float(abs(y[strata==s]-p[strata==s]).mean()),'magnitude_Spearman':correlation(y[strata==s],p[strata==s])} for s in np.unique(strata)}}
    if comparison is not None:result.update(relative_MAE_gain_vs_sequence=1-error/macro_error(y,comparison,strata),gain_vs_sequence_95CI=bounds(compare))
    return result

def model_gate(m,protocol):
    return bool(m['relative_MAE_gain_vs_line_median']>=protocol['model_required_relative_MAE_gain_vs_cell_and_line_medians'] and m['relative_MAE_gain_vs_cell_median']>=protocol['model_required_relative_MAE_gain_vs_cell_and_line_medians']
        and m['gain_vs_line_95CI'] is not None and m['gain_vs_line_95CI'][0]>0 and m['gain_vs_cell_95CI'] is not None and m['gain_vs_cell_95CI'][0]>0
        and m['macro_within_cell_line_Spearman'] is not None and m['macro_within_cell_line_Spearman']>=protocol['magnitude_Spearman_minimum'] and m['macro_within_line_Spearman_95CI'] is not None and m['macro_within_line_Spearman_95CI'][0]>0)

def pair_metrics(y,p,base,pairs,mask,groups,strata,raw,calibration,repetitions=2000,seed=20261108):
    pairs=pairs[np.all(mask[pairs],axis=1)];a,b=pairs[:,0],pairs[:,1]
    if not len(pairs):return {'pairs':0}
    dy=y[a]-y[b];dp=p[a]-p[b];db=base[a]-base[b];nonzero=dy!=0;error=float(abs(dy-dp).mean());zero=float(abs(dy).mean());be=float(abs(dy-db).mean());width=np.maximum(calibration['factor']*(raw[a]+raw[b]),calibration['absolute_error_floor']);keys=np.array([str(strata[i])+'|'+str(strata[j]) for i,j in pairs]);within=[]
    for key in np.unique(keys):
        ix=keys==key
        if ix.sum()>=30:
            r=correlation(dy[ix],dp[ix])
            if r is not None:within.append(r)
    result={'pairs':len(pairs),'genomic_1Mb_windows':len(np.unique(groups[a])),'difference_MAE':error,'zero_difference_MAE':zero,'training_median_difference_MAE':be,
            'relative_gain_vs_training_median_difference':1-error/be if be else None,'difference_Spearman':correlation(dy,dp),
            'macro_within_fixed_cell_pair_Spearman':float(np.mean(within)) if within else None,'fixed_cell_pair_strata':len(within),
            'ordering_accuracy_nonzero_measured_differences':float(np.mean(np.sign(dy[nonzero])==np.sign(dp[nonzero]))) if nonzero.any() else None,
            'direct_difference_interval_coverage':float(np.mean(abs(dy-dp)<=width)),'median_difference_interval_halfwidth':float(np.median(width))}
    windows=np.unique(groups[a]);result['difference_gain_vs_training_median_95CI']=None;result['difference_Spearman_95CI']=None;result['macro_within_fixed_cell_pair_Spearman_95CI']=None;result['ordering_accuracy_95CI']=None
    if len(windows)>=2:
        parts=[np.flatnonzero(groups[a]==g) for g in windows];rng=np.random.default_rng(seed);gain=[];rhos=[];within_rhos=[];orders=[]
        for iteration in range(repetitions):
            ix=np.concatenate([parts[j] for j in rng.integers(len(parts),size=len(parts))]);e=float(abs(dy[ix]-dp[ix]).mean());b=float(abs(dy[ix]-db[ix]).mean());gain.append(1-e/b if b else None);rhos.append(correlation(dy[ix],dp[ix]));nz=dy[ix]!=0;orders.append(float(np.mean(np.sign(dy[ix][nz])==np.sign(dp[ix][nz]))) if nz.any() else None);values=[]
            for key in np.unique(keys[ix]):
                selected=ix[keys[ix]==key]
                if len(selected)>=30:
                    r=correlation(dy[selected],dp[selected])
                    if r is not None:values.append(r)
            within_rhos.append(float(np.mean(values)) if values else None)
        result.update(difference_gain_vs_training_median_95CI=bounds(gain),difference_Spearman_95CI=bounds(rhos),macro_within_fixed_cell_pair_Spearman_95CI=bounds(within_rhos),ordering_accuracy_95CI=bounds(orders))
    return result

def main():
    start=time.monotonic();assert_previous_sources();lock_path=FINAL/'before_reserved_outcomes_lock.json';lock=json.loads(lock_path.read_text());assert lock['all_reserved_predictions_and_confidence_fixed'] and not lock['reserved_outcome_join_or_performance_calculated'];assert not (FINAL/'reserved_outcome_evaluation.json').exists()
    for path,h in lock['files'].items():assert sha(path)==h,path
    rows=list(csv_rows(DEST/'reserved_observation_metadata.csv.gz'));physical=list(csv_rows(DEST/'reserved_physical_sequence_bank.csv.gz'));lookup={}
    # The first numeric target join happens here, after both predictions were locked.
    with (SCREEN/'Tardaguila2026_TableS3.tsv').open() as handle:
        for r in csv.DictReader(handle,delimiter='\t'):
            key=('Tardaguila2026',r['SNP'],r['comparison'],r['window'],r['Cell_Type']);assert key not in lookup;lookup[key]=(float(r['log2Skew']),float(r['Skew_SE']))
    source_receipt=json.loads((SCREEN/'Siraj_unused_libraries/targeted_retrieval_receipt.json').read_text())
    for item in source_receipt['files']:
        with gzip.open(ROOT/item['path'],'rt') as handle:
            for r in csv.DictReader(handle,delimiter='\t'):
                if r['comb']!=r['SNP']+'_center_fwd_ref':continue
                key=('Siraj2026',item['library'],'chr'+r['SNP'].removeprefix('chr'),item['cell']);assert key not in lookup
                try:value=float(r['Log2Skew']);se=float(r['Skew_SE'])
                except ValueError:value=se=float('nan')
                lookup[key]=(value,se)
    y=[];ses=[];groups=[]
    for r in rows:
        key=('Tardaguila2026',r['source_SNP'],r['comparison'],r['window'],r['source_cell']) if r['source']=='Tardaguila2026' else ('Siraj2026',r['library'],r['original_variant_GRCh37'],r['source_cell'])
        assert key in lookup,key;value,se=lookup[key];y.append(abs(value));ses.append(se);s=physical[int(r['bank_index'])];groups.append(s['chromosome']+':'+str((int(s['position_hg38'])-1)//1000000))
    y=np.asarray(y,dtype='float64');ses=np.asarray(ses);groups=np.asarray(groups);strata=np.array([r['cell_line'] for r in rows]);valid=np.isfinite(y);source=np.array([r['source'] for r in rows]);library=np.array([r['library'] for r in rows]);known=np.array([r['exact_line_known']=='True' for r in rows]);assert valid.all(),'No silently dropped unavailable outcomes; explicit audit amendment needed.'
    settings=json.loads((FINAL/'final_settings_lock.json').read_text());protocol=settings['new_test_protocol'];confidence_policy=json.loads((FINAL/'confidence_policy_lock.json').read_text());masks={'blood_K562':(source=='Tardaguila2026')&(strata=='K562'),'OL33_three_known_lines':(source=='Siraj2026')&(library=='OL33')}
    predictions={};results={};pair_results={};confidence={}
    with np.load(FINAL/'reserved_measured_pair_identities.npz') as z:pairs={k:z[k] for k in z.files}
    for arm in ['sequence','AVI']:
        with np.load(FINAL/(arm+'_reserved_predictions.npz')) as z:p={k:z[k] for k in z.files}
        assert np.array_equal(p['row_index'],np.arange(len(rows)));predictions[arm]=p;results[arm]={};pair_results[arm]={};confidence[arm]={};package=joblib.load(FINAL/(arm+'_confidence.joblib'))
        for name,mask in masks.items():
            ix=np.flatnonzero(mask);comparison=predictions['sequence']['magnitude'][ix] if arm=='AVI' else None
            metrics=primary_metrics(y[ix],p['magnitude'][ix],p['baseline'][ix],p['cell_baseline'][ix],groups[ix],strata[ix],protocol['bootstrap_replicates'],protocol['bootstrap_seed'],comparison);metrics['all_accuracy_gates_pass']=model_gate(metrics,protocol);results[arm][name]=metrics
            pair_results[arm][name]={key:pair_metrics(y,p['magnitude'],p['baseline'],q,mask,groups,strata,p['raw_error_scale'],package['paired_difference_scale'][key],protocol['bootstrap_replicates'],protocol['bootstrap_seed']) for key,q in pairs.items()}
            accepted=mask&p['pretest_policy_accepted'];coverage=float(np.mean(abs(y[mask]-p['magnitude'][mask])<=p['interval_halfwidth'][mask]));c={'rows':int(mask.sum()),'nominal_target':.9,'observed_coverage':coverage,'median_halfwidth':float(np.median(p['interval_halfwidth'][mask])),'accepted_rows':int(accepted.sum()),'accepted_windows':len(np.unique(groups[accepted])),'accepted_fraction':float(accepted.sum()/mask.sum()),'useful_confident_answers_pass':False}
            if accepted.any():
                c['accepted_coverage']=float(np.mean(abs(y[accepted]-p['magnitude'][accepted])<=p['interval_halfwidth'][accepted]));c['accepted_MAE']=float(abs(y[accepted]-p['magnitude'][accepted]).mean());c['accepted_baseline_MAE']=float(abs(y[accepted]-p['baseline'][accepted]).mean())
            if accepted.sum()>=protocol['minimum_accepted_examples'] and c['accepted_windows']>=protocol['minimum_accepted_1Mb_windows']:
                ai=np.flatnonzero(accepted);c['accepted_accuracy']=primary_metrics(y[ai],p['magnitude'][ai],p['baseline'][ai],p['cell_baseline'][ai],groups[ai],strata[ai],protocol['bootstrap_replicates'],protocol['bootstrap_seed'])
                c['useful_confident_answers_pass']=bool(metrics['all_accuracy_gates_pass'] and coverage>=confidence_policy['required_fresh_observed_coverage'] and c['accepted_coverage']>=confidence_policy['required_fresh_observed_coverage'] and model_gate(c['accepted_accuracy'],protocol))
            confidence[arm][name]=c;print('RESERVED TEST EVALUATED',arm,name,'MAE',metrics['macro_MAE'],'gain',metrics['relative_MAE_gain_vs_line_median'],'coverage',coverage,flush=True)
    local_s=pair_results['sequence']['OL33_three_known_lines']['local'];local_a=pair_results['AVI']['OL33_three_known_lines']['local'];enough=local_s.get('pairs',0)>=protocol['nearby_minimum_pairs'] and local_s.get('genomic_1Mb_windows',0)>=protocol['nearby_minimum_1Mb_windows'];ranking=enough and all(local_a.get(k) is not None and local_s.get(k) is not None and local_a[k]>=local_s[k]-protocol['nearby_ordering_and_difference_Spearman_max_drop'] for k in ['ordering_accuracy_nonzero_measured_differences','difference_Spearman'])
    avi_gain=all(m['relative_MAE_gain_vs_sequence']>=protocol['AVI_required_relative_MAE_gain_vs_sequence'] and m['gain_vs_sequence_95CI'] is not None and m['gain_vs_sequence_95CI'][0]>0 for m in results['AVI'].values());avi_quality=all(m['all_accuracy_gates_pass'] for m in results['AVI'].values());retain=bool(avi_gain and avi_quality and ranking)
    unknown=(source=='Tardaguila2026')&~known;diagnostics={arm:{'unfamiliar_blood_lines':{'rows':int(unknown.sum()),'raw_MAE':float(abs(y[unknown]-p['magnitude'][unknown]).mean()),'pretest_accepted':int(p['pretest_policy_accepted'][unknown].sum()),'authorized_answers':0,'scope':'Raw diagnostics only. Unknown cell lines share broad blood context; no validated line-specific predictions.'},
         'OL27_reserved':{'rows':int(((source=='Siraj2026')&(library=='OL27')).sum()),'scope':'Reserved connected components from training-library source; not independent publication.'}} for arm,p in predictions.items()}
    joined=[]
    for i,r in enumerate(rows):joined.append({**r,'absolute_original_Log2Skew':y[i],'original_Skew_SE':ses[i],'genomic_1Mb_window':groups[i],
         **{arm+'_'+k:p[k][i].item() for arm,p in predictions.items() for k in ['magnitude','baseline','interval_halfwidth','pretest_policy_accepted']},'authorized_project_answer':False})
    write_gzip_csv(FINAL/'reserved_joined_measurements_and_predictions.csv.gz',joined)
    result={'primary_accuracy':results,'paired_difference_results':pair_results,'confidence':confidence,'other_diagnostics':diagnostics,
            'AVI_decision':{'retain_AVI':retain,'both_primary_relative_gain_and_CI_gates_pass':avi_gain,'both_primary_absolute_accuracy_gates_pass':avi_quality,'OL33_local_confirmation_sufficient':bool(enough),'OL33_local_ranking_preserved':bool(ranking),'scope':'Permissive composite AVI only; no gain can be assigned uniquely to AlphaGenome.'},
            'before_outcome_lock_sha256':sha(lock_path),'all_lock_hashes_verified_before_join':True,'source_effect':'Absolute original tile-level paired GLM Log2Skew; not detected/not-detected labels or model-generated truth.',
            'independent_publication_validation':'Tardaguila2026 only. OL33 is a held library within the existing Siraj publication/lab; 236 blood designs are 49 distinct variants, not 236 independent variants.',
            'AVI_composite_training_exposure_fully_excluded':False,
            'AVI_training_exposure_scope':'Public composite used as a fixed feature. The three author DNA readers have separately audited literal/pretraining-library exclusions; this does not prove absence from every AlphaGenome/AVI training or tuning resource.',
            'these_test_outcomes_now_revealed':True,'cannot_reuse_as_untouched_for_changed_model':True,'independent_genetic_score_or_human_disease_validation':False,'causal_variant_gene_cell_proof':False,
            'authorized_project_answers':0,'scientific_objective_complete':False,'Nature_Communications_genetic_manuscript_ready':False,'elapsed_seconds':time.monotonic()-start}
    write_json(FINAL/'reserved_outcome_evaluation.json',result);assert_previous_sources();print('RESERVED OUTCOME EVALUATION COMPLETE',json.dumps(result['AVI_decision']),flush=True)

if __name__=='__main__':main()

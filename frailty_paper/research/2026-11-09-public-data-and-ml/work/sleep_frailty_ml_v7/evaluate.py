"""One locked whole-study test. Both input arms and adverse results are reported."""
import gzip,argparse,subprocess,sys
from train import *
from scipy.stats import spearmanr
def correlations(y,p,groups):
    values=[]
    for g in np.unique(groups):
        ix=groups==g
        if ix.sum()>=5 and np.ptp(y[ix])>0:values.append(float(spearmanr(y[ix],p[ix]).statistic) if np.ptp(p[ix]) else 0.)
    return {'global_Spearman':float(spearmanr(y,p).statistic) if np.ptp(p) else 0.,'local_mean_Spearman':float(np.mean(values)) if values else 0.,'local_groups':len(values)}
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--wait-for-lock',action='store_true');args=parser.parse_args()
    while args.wait_for_lock and not (OUT/'final_model_and_prediction_lock.json').exists():time.sleep(10)
    if (OUT/'reserved_study_results.json').exists():
        if (OUT/'additional_neural_validation_protocol.json').exists() and not (OUT/'combined_model_recommendation.json').exists():subprocess.run([sys.executable,str(WORK/'evaluate_additional_neural.py')],env=environment(),check=True)
        print('ALREADY EVALUATED; NO REFIT OR RETEST',flush=True);return
    lockpath=OUT/'final_model_and_prediction_lock.json';lock=json.loads(lockpath.read_text());assert lock['all_predictions_locked_before_test_performance'];assert sha(DATA/'reserved_original_effects.csv.gz')==lock['reserved_effect_file_sha256'];assert_previous_intact()
    for arm in ['sequence_only','AVI_added']:assert sha(lock[arm]['blind_predictions']['file'])==lock[arm]['blind_predictions']['sha256']
    if (OUT/'additional_neural_validation_protocol.json').exists():subprocess.run([sys.executable,str(WORK/'freeze_extra_validation.py')],env=environment(),check=True)
    write_json(OUT/'reserved_evaluation_started.json',{'time':time.time(),'lock_sha256':sha(lockpath),'prediction_and_model_parameters_fixed':True,'no_test_dependent_refit_permitted':True})
    with gzip.open(DATA/'reserved_original_effects.csv.gz','rt') as f:effects=list(csv.DictReader(f))
    design=list(csv.DictReader((DATA/'reserved_physical_design.csv').open()));assert len(effects)==9602;d=load();fit=np.flatnonzero((d['split']<2)&(d['unit']==0));medians={c:float(np.median(d['effect'][fit[d['cell'][fit]==i]])) for i,c in enumerate(d['cell_names'])}
    # Consecutive2kb components for local ranking. One-Mb genomic windows for
    # paired uncertainty; these are not ancestry-specific LD blocks.
    loci={};last=None;g=-1
    for r in sorted(design,key=lambda r:(r['chromosome'],int(r['position_hg38']))):
        c=r['chromosome'];p=int(r['position_hg38'])
        if last is None or c!=last[0] or p-last[1]>2000:g+=1
        loci[int(r['design_index'])]=g;last=(c,p)
    predictions={}
    for arm in ['sequence_only','AVI_added']:
        predictions[arm]={(int(r['design_index']),r['cell']):r for r in csv.DictReader(open(lock[arm]['blind_predictions']['file']))}
    joined=[]
    for e in effects:
        i=int(e['design_index']);cell='HNP' if e['cell_group']=='neural_progenitor' else 'HEK293';row={'design_index':i,'variant_GRCh38':design[i]['variant_GRCh38'],'cell':cell,'cell_group':e['cell_group'],'measured_absolute_log2_effect':float(e['absolute_log2_effect']),'original_logFC':float(e['original_logFC']),'original_p':float(e['original_p']),'original_q':float(e['original_q']),'baseline_cell_median':medians[e['cell_group']],'local_group':loci[i],'genomic_window':design[i]['chromosome']+':'+str((int(design[i]['position_hg38'])-1)//1000000)}
        for arm in predictions:
            r=predictions[arm][(i,cell)]
            for field in ['predicted_absolute_log2_effect','lower90','upper90','rankable','max_abs_standardized_feature']:row[arm+'_'+field]=float(r[field])
            row[arm+'_unranked_reason']=r['unranked_reason'];row[arm+'_AVI_available']=r['AVI_available']
        joined.append(row)
    path=OUT/'reserved_all_observations.csv'
    with path.open('w') as f:w=csv.DictWriter(f,fieldnames=list(joined[0]));w.writeheader();w.writerows(joined)
    y=np.array([r['measured_absolute_log2_effect'] for r in joined]);cell=np.array([r['cell'] for r in joined]);groups=np.array([r['local_group'] for r in joined]);baseline=np.array([r['baseline_cell_median'] for r in joined]);avi_known=np.array([str(r['AVI_added_AVI_available'])=='1' for r in joined]);arms={};errors={};localscore={}
    for arm in predictions:
        p=np.array([r[arm+'_predicted_absolute_log2_effect'] for r in joined]);lo=np.array([r[arm+'_lower90'] for r in joined]);hi=np.array([r[arm+'_upper90'] for r in joined]);rankable=np.array([r[arm+'_rankable'] for r in joined],dtype=bool);errors[arm]=abs(y-p);bycell={}
        for c in np.unique(cell):
            ix=cell==c;err=abs(y[ix]-p[ix]);bm=float(np.mean(abs(y[ix]-baseline[ix])));accepted=ix&rankable
            subsets={}
            for label,mask in [('exact_AVI_available',avi_known),('AVI_missing',~avi_known)]:
                use=ix&mask;subsets[label]={'n':int(use.sum()),'MAE':float(np.mean(abs(y[use]-p[use]))) if use.any() else None,'Spearman':float(spearmanr(y[use],p[use]).statistic) if use.sum()>=3 and np.ptp(y[use]) and np.ptp(p[use]) else None}
            bycell[c]={'n':int(ix.sum()),'MAE':float(err.mean()),'RMSE':float(np.sqrt(np.mean((y[ix]-p[ix])**2))),'baseline_MAE':bm,'MAE_skill':1-float(err.mean())/bm,**correlations(y[ix],p[ix],groups[ix]),'interval_coverage':float(np.mean((y[ix]>=lo[ix])&(y[ix]<=hi[ix]))),'median_interval_width':float(np.median(hi[ix]-lo[ix])),'ranked_fraction':float(np.mean(rankable[ix])),'ranked_observations':int(accepted.sum()),'ranked_MAE':float(np.mean(abs(y[accepted]-p[accepted]))) if accepted.any() else None,'ranked_interval_coverage':float(np.mean((y[accepted]>=lo[accepted])&(y[accepted]<=hi[accepted]))) if accepted.any() else None,'all_abstained':bool(not accepted.any()),'identical_row_AVI_availability_subsets':subsets}
        common=sorted(set(int(r['design_index']) for r in joined if r['cell']=='HNP')&set(int(r['design_index']) for r in joined if r['cell']=='HEK293'));lookup={(r['design_index'],r['cell']):r for r in joined};actual=np.array([lookup[(i,'HNP')]['measured_absolute_log2_effect']-lookup[(i,'HEK293')]['measured_absolute_log2_effect'] for i in common]);estimated=np.array([lookup[(i,'HNP')][arm+'_predicted_absolute_log2_effect']-lookup[(i,'HEK293')][arm+'_predicted_absolute_log2_effect'] for i in common]);reference=np.array([lookup[(i,'HNP')]['baseline_cell_median']-lookup[(i,'HEK293')]['baseline_cell_median'] for i in common]);mae=float(np.mean(abs(actual-estimated)));bm=float(np.mean(abs(actual-reference)));arms[arm]={'cells':bycell,'macro_MAE':float(np.mean([r['MAE'] for r in bycell.values()])),'local_mean_Spearman':float(np.mean([r['local_mean_Spearman'] for r in bycell.values()])),'cell_difference':{'matched_variants':len(common),'MAE':mae,'baseline_MAE':bm,'MAE_skill':1-mae/bm,'Spearman':float(spearmanr(actual,estimated).statistic) if np.ptp(estimated) else 0.,'ordering_accuracy_including_near_ties':float(np.mean(np.where((estimated==0)|(actual==0),.5,np.sign(estimated)==np.sign(actual))))}}
    windows=np.array([r['genomic_window'] for r in joined]);unique=np.unique(windows);sums=np.zeros((len(unique),2,2));counts=np.zeros((len(unique),2));baseline_sums=np.zeros((len(unique),2))
    for k,w in enumerate(unique):
        for ci,c in enumerate(['HNP','HEK293']):
            ix=(windows==w)&(cell==c);counts[k,ci]=ix.sum()
            baseline_sums[k,ci]=np.sum(abs(y[ix]-baseline[ix]))
            for ai,arm in enumerate(['sequence_only','AVI_added']):sums[k,ci,ai]=errors[arm][ix].sum()
    rng=np.random.default_rng(284901);improvements=[];sampled_MAEs=[];baseline_improvements=[]
    for _ in range(2000):
        sampled=rng.integers(len(unique),size=len(unique));n=counts[sampled].sum(0)
        if np.any(n==0):continue
        means=sums[sampled].sum(0)/n[:,None];sampled_MAEs.append(means);bm=baseline_sums[sampled].sum(0)/n;baseline_improvements.append(bm[:,None]-means);improvements.append(float(np.mean(means[:,0]-means[:,1])))
    sampled_MAEs=np.array(sampled_MAEs)
    for ai,arm in enumerate(['sequence_only','AVI_added']):
        for ci,c in enumerate(['HNP','HEK293']):arms[arm]['cells'][c]['MAE_genomic_window_CI95']=np.quantile(sampled_MAEs[:,ci,ai],[.025,.975]).tolist();arms[arm]['cells'][c]['MAE_improvement_over_cell_median_CI95']=np.quantile(np.array(baseline_improvements)[:,ci,ai],[.025,.975]).tolist();arms[arm]['cells'][c]['genomic_1Mb_windows']=int(len(np.unique(windows[cell==c])))
    ci=np.quantile(improvements,[.025,.975]).tolist();relative=1-arms['AVI_added']['macro_MAE']/arms['sequence_only']['macro_MAE'];criteria={'at_least_one_percent_relative_macro_MAE_improvement':relative>=.01,'paired_genomic_window_bootstrap_improvement_CI_above_zero':ci[0]>0,'neither_cell_MAE_worsens_over_five_percent':all(arms['AVI_added']['cells'][c]['MAE']<=1.05*arms['sequence_only']['cells'][c]['MAE'] for c in ['HNP','HEK293']),'local_ranking_drop_at_most_point_zero_one':arms['AVI_added']['local_mean_Spearman']>=arms['sequence_only']['local_mean_Spearman']-.01};retain=all(criteria.values());result={'study':'McAfee2023,doi10.1016/j.xgen.2023.100404','observations':len(joined),'physical_designs':len(design),'arms':arms,'AVI_comparison':{'relative_macro_MAE_improvement':relative,'paired_MAE_improvement_CI95':ci,'bootstrap_replicates':len(improvements),'one_Mb_genomic_windows':len(unique),'criteria':criteria,'retain_AVI_added_input':retain,'recommended_input_arm':'AVI_added' if retain else 'sequence_only','scope':'Prespecified conditional input recommendation from this one study, not a new confirmatory test of the selected arm. AVI composite benefit cannot be attributed uniquely to AlphaGenome.'},'source_limitation':'Original published fitted magnitudes match current summaries; raw barcode QC/model-fitting not independently reconstructed. HEK293 differs from development HEK293T; HNP differs from H1-derived development NPCs. No participant outcome or causal validation.','uncertainty_limitation':'Study shifts are not exchangeable by assumption. Report actual outside coverage and abstention. Genomic1Mb window intervals do not establish ancestry-matched LD independence.','joined_observations':str(path),'joined_sha256':sha(path),'lock_sha256':sha(lockpath),'no_test_dependent_refitting':True};write_json(OUT/'reserved_study_results.json',result);print('LOCKED STUDY RESULTS',json.dumps(result['AVI_comparison']),flush=True)
    if (OUT/'additional_neural_validation_protocol.json').exists():subprocess.run([sys.executable,str(WORK/'evaluate_additional_neural.py')],env=environment(),check=True)
if __name__=='__main__':main()

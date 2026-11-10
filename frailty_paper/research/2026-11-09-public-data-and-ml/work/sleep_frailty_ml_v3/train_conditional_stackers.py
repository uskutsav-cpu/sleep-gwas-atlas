"""Finite development-only combination experiment, no reserved inference.

These are selection-stage development comparisons, not unbiased nested CV:
base checkpoints used their folds for stopping, and other-fold base fits can
have seen the current meta-validation fold. The untouched reserved test is
therefore essential. This limitation must not be hidden by the term OOF.
"""
import pickle,shutil,subprocess
from common import *
from conditional_stacking import fit_models,predict_models
FAMILIES={'stacked_contrast_no_avi':('contrast_no_avi',True),'stacked_contrast_avi_teacher':('contrast_avi_teacher',True),'conditional_avi_only':('baseline_full_avi_direct',False)}
def main():
    assert_no_reserved_evaluation();spec=OUT/'conditional_stacking_prespecification.json';assert spec.exists();path=OUT/'development_comparison.json';base_copy=OUT/'development_comparison_before_stacking.json'
    if not base_copy.exists():shutil.copy2(path,base_copy)
    base=json.loads(base_copy.read_text());assert len(base['aggregate'])==9;data,_,manifest,_,_,_=load();active=np.array([int(r['active_screen']) for r in manifest]);y=data['label'];study=data['study'];development=np.flatnonzero(data['chromosome']>=4);assert len(set(development))==len(development);lookup={r['variant_GRCh38']:r for r in csv.DictReader((OUT/'avi_exact_allele_scores.csv').open())};avi=np.array([[float(lookup[r['variant_GRCh38']]['avi_raw']),float(lookup[r['variant_GRCh38']]['avi_phred']),1.] if lookup[r['variant_GRCh38']]['status']=='exact_allele_match' else [np.nan,np.nan,0.] for r in manifest],dtype='float32');fold_results=[];new_records=[];final_artifacts={};predictions={};fold_ids=np.full(len(y),-1,dtype='int32')
    folds=[list(range(4,10)),list(range(10,16)),list(range(16,23))];R_checks=[]
    for k,val in enumerate(folds):fold_ids[np.isin(data['chromosome'],val)]=k
    for family,(base_family,include_dna) in FAMILIES.items():
        expert=np.full(len(y),np.nan);overall=expert.copy();base_names={}
        for fold in base['fold_results']:
            if fold['family']!=base_family:continue
            val=fold['fold'];k=folds.index(val);name=fold['record_name'];prefix=name.rsplit('_ensemble_',1)[0] if include_dna else 'development_baselines_val'+'-'.join(map(str,val));cp=OUT/(prefix+('_conditional_validation_predictions.npz' if include_dna else '_validation_predictions.npz'))
            with np.load(cp) as z:ix=z['row_index'];assert set(data['chromosome'][ix])==set(val);expert[ix]=z[name]
            with np.load(OUT/(prefix+'_validation_predictions.npz')) as z:assert np.array_equal(z['row_index'],ix);overall[ix]=z[name]
            base_names[str(k)]=name
        assert np.isfinite(expert[development]).all() and np.isfinite(overall[development]).all();all_conditional=np.full(len(y),np.nan)
        for k,val in enumerate(folds):
            tr=np.flatnonzero((fold_ids>=0)&(fold_ids!=k)&(active==1));va=np.flatnonzero(fold_ids==k);models=fit_models(expert,avi,y,study,tr,include_dna);pred=predict_models(models,expert[va],avi[va],study[va],include_dna);all_conditional[va]=pred;conditional={}
            for s,study_name in enumerate(STUDIES):
                use=(study[va]==s)&(active[va]==1);conditional[study_name]=metrics(y[va[use]],pred[use])
                inp=OUT/(family+'_R_val'+'-'.join(map(str,val))+f'_{s}.csv');target=inp.with_name(inp.stem+'_metrics.csv')
                with inp.open('w',newline='') as f:w=csv.writer(f);w.writerow(['row_id','label','score']);w.writerows(zip(va[use],y[va[use]].astype('int32'),pred[use]))
                subprocess.run(['Rscript',str(ROOT/'work/sleep_frailty_ml_v2/independent_metrics.R'),str(inp),str(target)],capture_output=True,text=True,check=True);check=list(csv.DictReader(target.open()))[0];errors={key:abs(float(check[key])-conditional[study_name][key]) for key in ['auroc','average_precision']};assert max(errors.values())<1e-10;R_checks.append({'family':family,'fold':val,'study':study_name,'errors':errors})
            allmetrics=study_metrics(y[va],overall[va],study[va]);name=family+'_val'+'-'.join(map(str,val));checkpoint=MODELS/(name+'.pkl')
            with checkpoint.open('wb') as f:pickle.dump(models,f)
            fr={'fold':val,'family':family,'record_name':name,'all_calls':allmetrics,'conditional_active':conditional,'macro_conditional_AP_skill':float(np.mean([m['AP_skill_over_prevalence'] for m in conditional.values()])),'variant_score_definition':'Learned conditional-active combiner of cross-fold DNA expert and AVI' if include_dna else 'AVI-only classifier fitted on conditional-active development rows','members':[],'artifact':str(checkpoint.relative_to(OUT)),'meta_training_n':len(tr),'base_family':base_family,'meta_validation_qualification':'Development optimization; base stopping and other-fold fitting prevent an unbiased nested-CV interpretation.'};fold_results.append(fr);np.savez_compressed(OUT/(name+'_validation_predictions.npz'),row_index=va,overall=overall[va],conditional=pred,label=y[va],study=study[va]);print('CONDITIONAL COMBINATION',family,val,'skill',fr['macro_conditional_AP_skill'],flush=True)
        # Final coefficients fit development-only cross-fold base scores. The
        # deployment base encoders are then refitted using the locked epoch rule.
        tr=development[active[development]==1];models=fit_models(expert,avi,y,study,tr,include_dna);checkpoint=MODELS/('conditional_stacking_development_full_'+family+'.pkl')
        with checkpoint.open('wb') as f:pickle.dump(models,f)
        final_artifacts[family]={'artifact':str(checkpoint.relative_to(OUT)),'sha256':sha(checkpoint),'base_family':base_family,'include_DNA_expert':include_dna,'training_rows':len(tr),'training_indices_sha256':hashlib.sha256(tr.tobytes()).hexdigest(),'base_validation_records':base_names,'study_fits':{STUDIES[s]:{'n':m['training_n'],'calls':m['training_calls'],'coefficients':m['classifier'].coef_.tolist(),'intercept':m['classifier'].intercept_.tolist()} for s,m in models.items()}}
        predictions[family]=all_conditional[development]
    aggregates=[]
    for family in FAMILIES:
        rr=[r for r in fold_results if r['family']==family];aggregates.append({'family':family,'completed_folds':3,'all_fold_mean_all_call_AP_skill':float(np.mean([r['all_calls']['macro_AP_skill'] for r in rr])),'all_fold_mean_conditional_AP_skill':float(np.mean([r['macro_conditional_AP_skill'] for r in rr])),'per_study_mean':{s:{'all_call_auroc':float(np.mean([r['all_calls']['studies'][s]['auroc'] for r in rr])),'all_call_AP':float(np.mean([r['all_calls']['studies'][s]['average_precision'] for r in rr])),'conditional_auroc':float(np.mean([r['conditional_active'][s]['auroc'] for r in rr])),'conditional_AP':float(np.mean([r['conditional_active'][s]['average_precision'] for r in rr]))} for s in STUDIES}})
    output={**base,'fold_results':base['fold_results']+fold_results,'aggregate':base['aggregate']+aggregates,'R_replays':base['R_replays']+R_checks,'base_comparison_sha256':sha(base_copy),'conditional_stacking_specification_sha256':sha(spec),'stacking_qualification':__doc__,'test_evaluated':False};path.write_text(json.dumps(output,indent=2)+'\n');(OUT/'conditional_stacking_final_artifacts.json').write_text(json.dumps({'families':final_artifacts,'specification_sha256':sha(spec),'development_comparison_sha256':sha(path),'training_code_sha256':sha(Path(__file__)),'reserved_test_evaluated':False},indent=2)+'\n');np.savez_compressed(OUT/'conditional_stacking_development_predictions.npz',row_index=development,**predictions);print('LEARNED CONDITIONAL COMBINATIONS COMPLETE; RESERVED TEST UNTOUCHED',flush=True)
if __name__=='__main__':main()

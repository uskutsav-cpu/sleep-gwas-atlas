"""Open the reserved outcomes only after models AND predictions are frozen."""
import zipfile,subprocess,datetime,collections,sys
from common import *
from metrics import point
from scipy.stats import spearmanr
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v3'))
from blocked_metrics import paired_block_bootstrap
def truth_bool(s):
    assert s in ['TRUE','FALSE'],s
    return s=='TRUE'
def finite(s):
    try:return np.isfinite(float(s))
    except (ValueError,TypeError):return False
def variant_class(s):
    # Exact hierarchy from the study's section4.R, with explicit unknowns.
    rsid=s.get('rsid','') not in ['', 'NA'];common=s.get('ref_dbsnp') in ['minor','major'] and s.get('alt_dbsnp') in ['minor','major']
    if rsid and common:return 'natural_GWAS_alleles'
    if rsid:return 'GWAS_position_other_alleles'
    if s.get('snp_perimeter')=='TRUE':return 'LD_region_other_SNVs'
    if s.get('class')=='valentina_individual':return 'individual_synthetic_candidates'
    if s.get('spaced_mut')=='TRUE':return 'synthetic_saturation'
    return 'unknown_annotation'
def main():
    lock_path=OUT/'final_model_lock.json';lock=json.loads(lock_path.read_text())
    for path,h in lock['files'].items():assert sha(path)==h,path
    for path,h in lock['model_checkpoints'].items():assert sha(path)==h,path
    start=OUT/'final_test_evaluation_started.json'
    if start.exists():assert json.loads(start.read_text())['lock_sha256']==sha(lock_path)
    else:write_json(start,{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'lock_sha256':sha(lock_path),'purpose':'First final-study outcome opening. All subsequent results are observed; no retraining/reselection or fresh-test reclassification.'})
    sealed=json.loads((SOURCES/'kosicki_sealed_outcome_manifest.json').read_text());entry=next(r for r in sealed if r['member']=='mpravista/data/input/final_tests.rds');source_dir=DATA/'reserved_outcome_source';source_dir.mkdir(exist_ok=True);raw=source_dir/'final_tests.rds'
    if not raw.exists():
        with zipfile.ZipFile(SOURCES/'kosicki_mpravista.zip') as z:raw.write_bytes(z.read(entry['member']))
    assert sha(raw)==entry['sha256'];csv_path=source_dir/'final_tests_atomic.csv';db=source_dir/'dbAnnot.rds'
    if not db.exists():
        with zipfile.ZipFile(SOURCES/'kosicki_mpravista.zip') as z:db.write_bytes(z.read('mpravista/data/raw/dbAnnot.rds'))
    if not csv_path.exists():subprocess.run(['Rscript',str(WORK/'extract_reserved_outcomes.R'),str(raw),str(csv_path),str(start),str(db)],check=True,env=r_environment())
    all_rows=list(csv.DictReader(csv_path.open()));names=[r['name'] for r in all_rows];assert len(names)==len(set(names));lookup={r['name']:r for r in all_rows};design=list(csv.DictReader((DATA/'kosicki_blind_design.csv').open()));excluded=collections.Counter();keep=[];measurements=[]
    for i,r in enumerate(design):
        if r['name'] not in lookup:excluded['no_source_outcome_record']+=1;continue
        s=lookup[r['name']]
        if s['keep']!='TRUE':excluded['variant_failed_source_QC']+=1;continue
        if s['label_class']!='snp':excluded['source_not_SNV_label_class']+=1;continue
        ref=lookup.get(r['reference_name'])
        if ref is None or ref['keep']!='TRUE':excluded['missing_or_failed_QC_reference_construct']+=1;continue
        if 'BC_count' in s and (not finite(s['BC_count']) or float(s['BC_count'])<15 or not finite(ref.get('BC_count')) or float(ref['BC_count'])<15):excluded['source_variant_or_reference_fewer_than15barcodes']+=1;continue
        if not all(finite(s[k]) for k in ['padj','pval','zscore_diff']) or not 0<=float(s['padj'])<=1 or not 0<=float(s['pval'])<=1:excluded['unassessable_variant_outcome']+=1;continue
        called=truth_bool(s['sign']);assert called==(float(s['padj'])<.05)
        if all(finite(s.get(k)) for k in ['zscore_neg','zscore_ref']):
            assert abs(float(s['zscore_diff'])-(float(s['zscore_neg'])-float(s['zscore_ref'])))<1e-10
            if finite(ref.get('zscore_neg')) and abs(float(s['zscore_ref'])-float(ref['zscore_neg']))>1e-10:excluded['reference_readout_alignment_disagreement']+=1;continue
        keep.append(i);measurements.append({'label':int(called),'absolute_measured_effect':abs(float(s['zscore_diff'])),'variant_class':variant_class(s),'source_adjusted_p':float(s['padj']),'source_p':float(s['pval']),'source_var_id':s.get('var_id',''),'source_class':s.get('class',''),'source_spaced_mut':s.get('spaced_mut',''),'source_rsid':s.get('rsid','')})
    ix=np.array(keep,dtype='int32');assert len(ix)>0;y=np.array([r['label'] for r in measurements],dtype='int32');effect=np.array([r['absolute_measured_effect'] for r in measurements]);retained=[design[i] for i in ix]
    with np.load(OUT/'kosicki_blinded_predictions.npz') as z:
        assert z['name'].tolist()==[r['name'] for r in design];models={k.removeprefix('model::'):z[k][ix] for k in z.files if k.startswith('model::')};confidence={k:z[k][ix] for k in ['calibrated_call_probability','set_includes0','set_includes1','status','ranked_score']}
    blocks=np.array([r['chromosome']+':'+str(int(r['position_hg38'])//1000000) for r in retained]);tiles=np.array([r['reference_name'] for r in retained]);results={}
    for name,p in models.items():
        local_rho=[];local_auc=[];local_ap=[]
        for tile in np.unique(tiles):
            use=tiles==tile
            if use.sum()>=10 and len(np.unique(effect[use]))>1 and len(np.unique(p[use]))>1:local_rho.append(float(spearmanr(effect[use],p[use]).statistic))
            if use.sum()>=10 and y[use].sum()>=2 and (use.sum()-y[use].sum())>=2:
                q=point(y[use],p[use]);local_auc.append(q['AUROC']);local_ap.append(q['AP'])
        results[name]={'all_pairs':point(y,p),'effect_magnitude_Spearman':float(spearmanr(effect,p).statistic),'within_reference_tile_mean_Spearman':float(np.mean(local_rho)) if local_rho else None,'within_reference_tile_Spearman_groups':len(local_rho),'within_reference_tile_mean_AUROC':float(np.mean(local_auc)) if local_auc else None,'within_reference_tile_mean_AP':float(np.mean(local_ap)) if local_ap else None,'within_reference_tile_binary_groups':len(local_auc)}
    classes=np.array([r['variant_class'] for r in measurements]);subgroups={}
    for label in np.unique(classes):
        use=classes==label;subgroups[str(label)]={'pairs':int(use.sum()),'calls':int(y[use].sum()),'unique_genomic_contrasts':len({retained[i]['variant_GRCh38'] for i in np.flatnonzero(use)}),'models':{name:point(y[use],p[use]) for name,p in models.items()},'uncertainty_note':'Exploratory subgroup precision is limited, especially when few effects are detected. These are not sleep/frailty GWAS replications.'}
    ranked=confidence['status']=='ranked';positive=ranked&confidence['set_includes1'];negative=ranked&confidence['set_includes0'];covered=np.where(y==1,confidence['set_includes1'],confidence['set_includes0']);correct=(positive&(y==1))|(negative&(y==0));prob=confidence['calibrated_call_probability'];uncertainty={'total':len(y),'ranked':int(ranked.sum()),'ranked_fraction':float(ranked.mean()),'unranked_reasons':dict(collections.Counter(confidence['status'].tolist())),'prediction_set_label_coverage':float(covered.mean()),'coverage_among_detected':float(covered[y==1].mean()),'coverage_among_non_detected':float(covered[y==0].mean()),'ranked_singleton_accuracy':float(correct[ranked].mean()) if ranked.any() else None,'ranked_predicted_detected':int(positive.sum()),'ranked_detected_precision':float(y[positive].mean()) if positive.any() else None,'ranked_detected_recall':float(np.sum(positive&(y==1))/y.sum()),'calibrated_probability_Brier':float(np.mean((prob-y)**2)),'caution':'Nominal90%prediction sets are not90%positive predictive probabilities, and source-shift coverage is empirical. Unranked rows have no ranked_score.'}
    with (OUT/'reserved_study_predictions.csv').open('w',newline='') as f:
        fields=['row_id','name','reference_name','variant_GRCh38','chromosome','position_hg38','label','absolute_measured_effect','variant_class','status','ranked_score','calibrated_laboratory_call_probability',*models];w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for j,(r,m) in enumerate(zip(retained,measurements)):w.writerow({'row_id':int(ix[j]),**{k:r[k] for k in fields[1:6]},'label':m['label'],'absolute_measured_effect':m['absolute_measured_effect'],'variant_class':m['variant_class'],'status':str(confidence['status'][j]),'ranked_score':float(confidence['ranked_score'][j]) if np.isfinite(confidence['ranked_score'][j]) else '', 'calibrated_laboratory_call_probability':float(prob[j]),**{k:float(p[j]) for k,p in models.items()}})
    np.savez_compressed(OUT/'reserved_analysis_arrays.npz',design_index=ix,label=y,effect=effect,blocks=blocks,tiles=tiles,variant_class=classes,**{'model::'+k:p for k,p in models.items()},**confidence)
    bootstrap=paired_block_bootstrap(y,models,blocks,repetitions=1000,seed=20261107);write_json(OUT/'reserved_block_bootstrap.json',bootstrap)
    src=OUT/'reserved_R_metric_input.csv';dst=OUT/'reserved_R_metrics.csv'
    with src.open('w',newline='') as f:w=csv.writer(f);w.writerow(['row_id','label',*models]);w.writerows(zip(ix,y,*models.values()))
    subprocess.run(['Rscript',str(ROOT/'work/sleep_frailty_ml_v2/independent_metrics.R'),str(src),str(dst)],check=True,capture_output=True,text=True,env=r_environment());rchecks=[]
    for r in csv.DictReader(dst.open()):
        q=results[r['model']]['all_pairs'];error=max(abs(float(r['auroc'])-q['AUROC']),abs(float(r['average_precision'])-q['AP']));assert error<1e-12;rchecks.append({'model':r['model'],'max_error':error})
    report={'source_rows':len(all_rows),'source_QC_passed_SNV_rows':sum(r['keep']=='TRUE' and r['label_class']=='snp' for r in all_rows),'verified_SNV_design_rows':len(design),'eligible_pairs':len(y),'unique_genomic_contrasts':len({r['variant_GRCh38'] for r in retained}),'reference_tiles':len(np.unique(tiles)),'detected_effect_pairs':int(y.sum()),'exclusions':dict(excluded),'models':results,'variant_class_subgroups':subgroups,'uncertainty':uncertainty,'independent_R_checks':rchecks,'model_lock_sha256':sha(lock_path),'outcome_export_sha256':sha(csv_path),'predictions_sha256':sha(OUT/'reserved_study_predictions.csv'),'frozen_prediction_sha256':sha(OUT/'kosicki_blinded_predictions.npz'),'scope':'New held-out study outcomes, evaluated after complete model/threshold/prediction lock. Distinct measurements and publication, but both main neuronal sources use WTC11-derivedNGN2neurons; donor independence is not established. Not patient insomnia/frailty prediction, genetic-score validation, causation or an unseen-reference-genome test.','next':'Report every prespecified result without retraining/reselection on this now-observed study.'}
    write_json(OUT/'reserved_study_results.json',report);print('FROZEN NEW-STUDY EVALUATION COMPLETE',len(y),'pairs;',int(y.sum()),'detected;',len(rchecks)*2,'independentRmetrics',flush=True)
if __name__=='__main__':main()

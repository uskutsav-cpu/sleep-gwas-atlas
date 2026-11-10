"""Versioned extra-study test; original frozen evaluator and its source failure preserved."""
import zipfile,subprocess,datetime,collections,sys,gzip
from common import *
from metrics import point
from scipy.stats import spearmanr
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v3'))
from blocked_metrics import paired_block_bootstrap
def finite(s):
    try:return np.isfinite(float(s))
    except (ValueError,TypeError):return False
def main():
    lock_path=OUT/'final_model_lock.json';lock=json.loads(lock_path.read_text())
    for path,h in lock['files'].items():assert sha(path)==h,path
    for path,h in lock['model_checkpoints'].items():assert sha(path)==h,path
    additional_lock=OUT/'additional_study_validation_lock.json';additional=json.loads(additional_lock.read_text())
    assert additional['original_model_lock_sha256']==sha(lock_path)
    for p,h in additional['files'].items():assert sha(p)==h,p
    start=OUT/'additional_test_evaluation_started.json'
    if start.exists():assert json.loads(start.read_text())['lock_sha256']==sha(lock_path)
    else:write_json(start,{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'lock_sha256':sha(lock_path),'additional_validation_lock_sha256':sha(additional_lock),'purpose':'First final-study outcome opening. All subsequent results are observed; no retraining/reselection or fresh-test reclassification.'})
    csv_path=SOURCES/'recovery_IGVFFI4642CBGA.gz'
    with gzip.open(csv_path,'rt') as f:all_rows=list(csv.DictReader(f,delimiter='\t'))
    counts=collections.Counter(r['variant_id'] for r in all_rows);lookup={r['variant_id']:r for r in all_rows if counts[r['variant_id']]==1}
    bed={};bed_counts=collections.Counter()
    with gzip.open(SOURCES/'recovery_IGVFFI7902MFBK.gz','rt') as f:
        for line in f:
            fields=line.rstrip().split('\t');assert len(fields)==19;bed_counts[fields[3]]+=1;bed[fields[3]]=fields
    design=list(csv.DictReader((DATA/'additional_study_blind_design.csv').open()));excluded=collections.Counter();keep=[];measurements=[]
    for i,r in enumerate(design):
        key=r['SPDI'];s=lookup.get(key);b=bed.get(key)
        if s is None or b is None or bed_counts[key]!=1:excluded['missing_or_duplicate_canonical_variant_outcome']+=1;continue
        if b[0]!=r['chromosome'] or int(b[1])+1!=int(r['position_hg38']) or int(b[2])!=int(b[1])+1:excluded['genomic_coordinate_disagreement']+=1;continue
        shared=[k for k in list(s) if k!='variant_id']
        if any(str(s[k])!=str(v) for k,v in zip(shared,b[6:])):excluded['current_genomic_BED_vs_variant_TSV_disagreement']+=1;continue
        if not all(finite(s[k]) for k in ['log2FoldChange','minusLog10PValue','minusLog10QValue']) or float(s['minusLog10PValue'])<0 or float(s['minusLog10QValue'])<0:excluded['missing_or_unassessable_measured_effect_p_q']+=1;continue
        if abs(float(s['minusLog10QValue'])+np.log10(.05))<=.00005:excluded['source_qvalue_rounding_interval_contains_0p05']+=1;continue
        ra,aa=key.split(':')[-2:]
        if s['refAllele']!=ra or s['altAllele']!=aa:excluded['canonical_forward_genomic_allele_disagreement']+=1;continue
        called=float(s['minusLog10QValue'])>-np.log10(.05)
        keep.append(i);measurements.append({'label':int(called),'absolute_measured_effect':abs(float(s['log2FoldChange'])),'variant_class':'CNN_selected_test_SNV','source_adjusted_p':float(10**(-float(s['minusLog10QValue']))),'source_p':float(10**(-float(s['minusLog10PValue']))),'source_var_id':key,'source_class':'test','source_spaced_mut':'','source_rsid':''})
    ix=np.array(keep,dtype='int32');assert len(ix)>0;y=np.array([r['label'] for r in measurements],dtype='int32');effect=np.array([r['absolute_measured_effect'] for r in measurements]);retained=[design[i] for i in ix]
    with np.load(OUT/'additional_study_blinded_predictions.npz') as z:
        assert z['name'].tolist()==[r['name'] for r in design];models={k.removeprefix('model::'):z[k][ix] for k in z.files if k.startswith('model::')};confidence={k:z[k][ix] for k in ['calibrated_call_probability','set_includes0','set_includes1','status','ranked_score']}
    blocks=np.array([r['chromosome']+':'+str(int(r['position_hg38'])//1000000) for r in retained]);tiles=np.array([r['reference_name'] for r in retained]);results={};tile_indices=collections.defaultdict(list)
    for j,t in enumerate(tiles):tile_indices[str(t)].append(j)
    for name,p in models.items():
        local_rho=[];local_auc=[];local_ap=[]
        for indices in tile_indices.values():
            use=np.array(indices,dtype='int32')
            if len(use)>=10 and len(np.unique(effect[use]))>1 and len(np.unique(p[use]))>1:local_rho.append(float(spearmanr(effect[use],p[use]).statistic))
            if len(use)>=10 and y[use].sum()>=2 and (len(use)-y[use].sum())>=2:
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
    report={'source_rows':len(all_rows),'study':'IGVFDS6149YGOL HepG2 / IGVFSM9009DVDG','canonical_current_variant_rows':len(lookup),'verified_SNV_design_rows':len(design),'eligible_pairs':len(y),'unique_genomic_contrasts':len({r['variant_GRCh38'] for r in retained}),'reference_tiles':len(np.unique(tiles)),'detected_effect_pairs':int(y.sum()),'exclusions':dict(excluded),'models':results,'variant_class_subgroups':subgroups,'uncertainty':uncertainty,'independent_R_checks':rchecks,'model_lock_sha256':sha(lock_path),'additional_validation_lock_sha256':sha(additional_lock),'outcome_export_sha256':sha(csv_path),'predictions_sha256':sha(OUT/'reserved_study_predictions.csv'),'frozen_prediction_sha256':sha(OUT/'additional_study_blinded_predictions.npz'),'scope':'Additional whole-assay HepG2 measurements, all models/thresholds fixed before original test opening and additional predictions frozen before any additional-study performance calculation. Entire study source distinct from training, but same established HepG2 cell line and overlapping laboratory network;3technical replicates do not establish independent donors. CNN-selected variants and source processing versionv0.0.7. No clinical disease, genetic-score or causal validation.','source_recovery_qualification':'Original Kosicki canonical binary test failed because all SNP statistical fields were missing. Actual continuous measurements were recovered separately; this extra complete source supplies a separate binary/abstention test. No model, threshold or primary selection changed.','next':'Report every prespecified result without retraining/reselection on this now-observed study.'}
    write_json(OUT/'reserved_study_results.json',report);print('FROZEN NEW-STUDY EVALUATION COMPLETE',len(y),'pairs;',int(y.sum()),'detected;',len(rchecks)*2,'independentRmetrics',flush=True)
if __name__=='__main__':main()

"""One explicitly locked reserved evaluation, with fixed paired uncertainty.

Only run after final selection, refits and verification. No test-driven redesign
or reselection is allowed. --resume can finish the same immutable evaluation.
"""
import argparse,datetime,subprocess,json,csv,numpy as np
from common import *
from inference_records import StoredFeatureInference
from blocked_metrics import paired_block_bootstrap
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--registry-sha256',required=True);parser.add_argument('--resume',action='store_true');args=parser.parse_args();path=OUT/'final_model_registry.json';registry=json.loads(path.read_text());lock=json.loads((OUT/'final_artifact_lock.json').read_text());verify=json.loads((OUT/'final_fit_verification.json').read_text());assert args.registry_sha256==sha(path)==lock['registry_sha256']==verify['registry_sha256'];assert lock['verification_sha256']==sha(OUT/'final_fit_verification.json') and verify['no_reserved_test_scores_evaluated'];assert not (OUT/'reserved_evaluation_report.json').exists(),'Do not repeat a completed evaluation'
    plan=json.loads((OUT/'final_training_locked_plan.json').read_text());assert sha(OUT/'final_training_locked_plan.json')==registry['plan_sha256']==lock['final_training_plan_sha256']==verify['final_plan_sha256'];assert plan['protocol_sha256']==sha(OUT/'final_selection_protocol.json');assert registry['source_code_sha256']==sha(WORK/'final_fit.py')
    exposure=json.loads((OUT/'reserved_exposure_audit.json').read_text());assert lock['reserved_exposure_audit_sha256']==sha(OUT/'reserved_exposure_audit.json')
    for stem in ['salomon','siraj']:assert exposure['source_manifest_sha256'][stem]==sha(WORK/(stem+'_manifest_local.csv'))
    started=OUT/'reserved_evaluation_started.json';signature={'registry_sha256':sha(path),'final_plan_sha256':sha(OUT/'final_training_locked_plan.json'),'selection_protocol_sha256':sha(OUT/'final_selection_protocol.json'),'reserved_exposure_audit_sha256':sha(OUT/'reserved_exposure_audit.json'),'evaluation_code_sha256':sha(Path(__file__)),'bootstrap_code_sha256':sha(WORK/'blocked_metrics.py')}
    if started.exists():
        assert args.resume,'Use --resume only to finish the existing immutable evaluation';prior=json.loads(started.read_text());assert all(prior[k]==v for k,v in signature.items())
    else:
        assert not args.resume
        with started.open('x') as f:json.dump({**signature,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'First reserved v3 model evaluation on Siraj fixed-library chromosomes1–3. Salomon study data are previously used development evidence. No further v3 fitting/selection permitted.'},f,indent=2)
    engine=StoredFeatureInference();data=engine.data;ix=np.flatnonzero(data['chromosome']<=3);assert set(data['chromosome'][ix])=={1,2,3};records=registry['records'];selected=sorted(set(registry['role_model_names'].values()))
    for r in records.values():
        if r['artifact']:assert sha(OUT/r['artifact'])==r['checkpoint_sha256']
    cache=OUT/'reserved_test_predictions.npz';cache_receipt=OUT/'reserved_test_prediction_receipt.json'
    if cache.exists():
        receipt=json.loads(cache_receipt.read_text());assert receipt['signature']==signature and receipt['prediction_sha256']==sha(cache)
        with np.load(cache) as z:assert np.array_equal(z['row_index'],ix);predictions={k:z[k] for k in z.files if k!='row_index'}
    else:
        predictions={}
        for name in selected:
            p=engine.predict(records[name],ix,records);predictions[name+'__overall']=p['overall'];predictions[name+'__active_variant']=p['conditional_scores'][:,1] if 'conditional_scores' in p else p['overall'];print('Reserved predictions computed for locked model',name,flush=True)
        tmp=cache.with_name(cache.stem+'.tmp.npz');np.savez_compressed(tmp,row_index=ix,**predictions);tmp.replace(cache);cache_receipt.write_text(json.dumps({'signature':signature,'prediction_sha256':sha(cache),'row_indices_sha256':hashlib.sha256(ix.tobytes()).hexdigest(),'selected_model_names':selected},indent=2)+'\n')
    assert set(predictions)=={n+s for n in selected for s in ['__overall','__active_variant']};active=engine.features.active[ix];coordinates=np.array([int(engine.manifest[i]['position_hg38']) for i in ix]);results={};R_checks=[];strict=np.isin(ix,exposure['retained_Siraj_global_row_indices']);assert int(strict.sum())==exposure['strict_Siraj_reserved_n']
    for s,study in enumerate(STUDIES):
        results[study]={'evidence_status':'Reserved fixed-library chromosome evaluation; foundation-genome exposure not excluded' if s==1 else 'Previously used study; not an independent fresh-study test','strata':{}}
        study_mask=data['study'][ix]==s
        full_comparison={}
        if s==1:
            for label,mask0,suffix0 in [('all_calls',study_mask,'__overall'),('measured_active_regions',study_mask&(active==1),'__active_variant')]:full_comparison[label]={n:metrics(data['label'][ix[mask0]],predictions[n+suffix0][mask0]) for n in selected}
            results[study]['original_full_reserved_set_comparison']={'n':int(study_mask.sum()),'additional_prior_project_region_exclusions':exposure['additional_excluded_n'],'point_metrics':full_comparison,'qualification':'Original full set includes regions close to previous project molecular studies. Strict region-filtered set is the primary reserved evaluation; eligibility was refined before any reserved prediction.'};study_mask=study_mask&strict
        for stratum,mask,suffix in [('all_calls',study_mask,'__overall'),('measured_active_regions',study_mask&(active==1),'__active_variant')]:
            y=data['label'][ix[mask]].astype('int32');pp={n:predictions[n+suffix][mask] for n in selected};metrics_by_model={n:metrics(y,p) for n,p in pp.items()};primary_blocks=data['chromosome'][ix[mask]]*10000+(coordinates[mask]-1)//1000000;sensitivity_blocks=data['chromosome'][ix[mask]]*10000+(coordinates[mask]-1)//5000000
            uncertainty=paired_block_bootstrap(y,pp,primary_blocks,repetitions=1000,seed=2311+s*100+(stratum!='all_calls'));sensitivity=paired_block_bootstrap(y,pp,sensitivity_blocks,repetitions=500,seed=7811+s*100+(stratum!='all_calls'))
            inp=OUT/(f'reserved_R_{s}_{stratum}.csv');target=inp.with_name(inp.stem+'_metrics.csv')
            with inp.open('w',newline='') as f:
                w=csv.writer(f);w.writerow(['row_id','label',*selected]);w.writerows(zip(ix[mask],y,*[pp[n] for n in selected]))
            subprocess.run(['Rscript',str(ROOT/'work/sleep_frailty_ml_v2/independent_metrics.R'),str(inp),str(target)],check=True,capture_output=True,text=True)
            for r in csv.DictReader(target.open()):
                name=r['model'];error={k:abs(float(r[k])-metrics_by_model[name][k]) for k in ['auroc','average_precision']};assert max(error.values())<1e-10;R_checks.append({'study':study,'stratum':stratum,'model':name,'errors':error})
            results[study]['strata'][stratum]={'point_metrics':metrics_by_model,'primary_1Mb_block_bootstrap':uncertainty,'sensitivity_5Mb_block_bootstrap':sensitivity,'active_subset_definition':'Original measured-active assay screen used only for evaluation strata; never supplied at prediction time' if stratum!='all_calls' else None};print('Finished reserved evaluation',study,stratum,'n',len(y),flush=True)
    output={'signature':signature,'role_model_names':registry['role_model_names'],'results':results,'independent_R_checks':R_checks,'prediction_sha256':sha(cache),'scope':'Predicting laboratory assay calls/variant effects; not patient insomnia/frailty outcomes, causal genes, genetic correlation, polygenic-score validation, independent donors or clinical generalization. No test-based model revision/reselection.','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};(OUT/'reserved_evaluation_report.json').write_text(json.dumps(output,indent=2)+'\n');print('LOCKED RESERVED EVALUATION COMPLETE; DO NOT RESELECT ON THESE RESULTS',flush=True)
if __name__=='__main__':main()

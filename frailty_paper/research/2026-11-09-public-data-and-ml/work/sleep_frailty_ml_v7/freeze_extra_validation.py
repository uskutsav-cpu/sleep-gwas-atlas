"""Lock every new neuronal prediction before opening either outside performance."""
from finalize import *
from sequence_features import Encoder
def main():
    path=OUT/'additional_neural_prediction_lock.json'
    if path.exists():
        lock=json.loads(path.read_text())
        for r in lock['arms'].values():assert sha(r['file'])==r['sha256']
        return
    assert not (OUT/'reserved_evaluation_started.json').exists();assert not (OUT/'additional_neural_results.json').exists();protocol=json.loads((OUT/'additional_neural_validation_protocol.json').read_text());assert sha(DATA/'neuroimmune_physical_design_verified.csv')==protocol['design_sha256'];assert sha(DATA/'neuroimmune_verified_candidate_effects.csv')==protocol['effects_sha256'];primary=json.loads((OUT/'final_model_and_prediction_lock.json').read_text());extra=json.loads((OUT/'extra_AVI_inference_receipt.json').read_text());assert sha(extra['output'])==extra['sha256']
    rows=[r for r in csv.DictReader((DATA/'neuroimmune_physical_design_verified.csv').open()) if not r['additional_test_overlap_exclusion']];assert len(rows)==protocol['physical_predictions_to_freeze'];x=Encoder().features([r['reference'] for r in rows],[r['alternate'] for r in rows],[int(r['changed_position']) for r in rows]);features_path=DATA/'additional_neural_feature_bank.npy';np.save(features_path,x);d=load();test=dict(d);n=len(rows);cell=d['cell_names'].index('induced_excitatory_neuron');test.update({'x':x,'n':n,'bank_index':np.arange(n,dtype='int32'),'unit':np.zeros(n,dtype='int8'),'cell':np.full(n,cell,dtype='int16'),'line':np.full(n,-1,dtype='int16')});lookup={r['variant_GRCh38']:r for r in csv.DictReader(open(extra['output']))};avi=[]
    for r in rows:
        value=lookup[r['variant_GRCh38']];known=int(value['available']);avi.append([float(value['raw']) if known else np.nan,float(value['phred']) if known else np.nan,known])
    avi=np.array(avi,dtype='float32');arms={}
    for arm in ['sequence_only','AVI_added']:
        records=primary[arm]['records'];cal=json.loads((OUT/(arm+'_calibration.json')).read_text());p=predict_ensemble(test,np.arange(n),records,avi if arm=='AVI_added' else None);tx=Path(records[0]['record']['transform']);t=transform_load(tx);maxz=np.max(abs((x-t['mean'])/t['sd']),axis=1);scale=cal['cell_error_scale'].get(str(cell),cal['global_error_scale']);half=scale*cal['worst_study_standardized_q90'];out=OUT/f'additional_neural_{arm}_blind_predictions.csv'
        with out.open('w') as f:
            writer=csv.DictWriter(f,fieldnames=['design_index','source_ID','variant_GRCh38','cell_group','baseline_context_projection','lower90','upper90','frozen_rule_rankable','unranked_reason','max_abs_standardized_feature','AVI_available']);writer.writeheader()
            for i,r in enumerate(rows):
                lo=max(0.,p[i]-half);hi=p[i]+half;reasons=[]
                if cal['cell_training_study_support'].get(str(cell),0)<2:reasons.append('cell_group_supported_by_fewer_than_two_training_studies')
                if maxz[i]>12:reasons.append('sequence_features_outside_frozen_domain_limit')
                if hi-lo>cal['maximum_interval_width']:reasons.append('interval_too_wide')
                writer.writerow({'design_index':r['design_index'],'source_ID':r['source_ID'],'variant_GRCh38':r['variant_GRCh38'],'cell_group':'induced_excitatory_neuron','baseline_context_projection':p[i],'lower90':lo,'upper90':hi,'frozen_rule_rankable':int(not reasons),'unranked_reason':';'.join(reasons),'max_abs_standardized_feature':maxz[i],'AVI_available':int(avi[i,2])})
        arms[arm]={'file':str(out),'sha256':sha(out),'rows':n}
    lock={'arms':arms,'feature_bank_sha256':sha(features_path),'extra_AVI_receipt_sha256':sha(OUT/'extra_AVI_inference_receipt.json'),'protocol_sha256':sha(OUT/'additional_neural_validation_protocol.json'),'primary_model_lock_sha256':sha(OUT/'final_model_and_prediction_lock.json'),'source_effect_file_sha256':protocol['effects_sha256'],'source_design_sha256':protocol['design_sha256'],'source_effect_values_read_in_this_prediction_script':False,'any_primary_or_additional_model_performance_seen':False,'stimulated_conditions_authorized_for_prediction':False};write_json(path,lock);print('ADDITIONAL NEURAL PREDICTIONS LOCKED BEFORE BOTH OUTSIDE TESTS',n,flush=True)
if __name__=='__main__':main()

"""Conservative research paired-DNA magnitude scorer, not a disease-risk model."""
import argparse,gzip
from finalize import *
from sequence_features import Encoder
def build_registry():
    d=load();lock=json.loads((OUT/'final_model_and_prediction_lock.json').read_text());result=json.loads((OUT/'reserved_study_results.json').read_text());combined=json.loads((OUT/'combined_model_recommendation.json').read_text());arm=combined['recommended_input_arm'];cal=json.loads((OUT/(arm+'_calibration.json')).read_text())
    registry={'recommended_arm':arm,'records':lock[arm]['records'],'calibration':cal,'cell_names':d['cell_names'],'cell_line_names':d['cell_line_names'],'family_by_cell':d['family_by_cell'].tolist(),'family_names':d['family_names'],'n_context_flags':d['n_context_flags'],'outside_study_results_sha256':sha(OUT/'reserved_study_results.json'),'combined_recommendation_sha256':sha(OUT/'combined_model_recommendation.json'),'cell_group_authorized_for_automated_research_ranking':combined['cell_group_authorized_for_automated_research_ranking'],'external_validation_by_cell_group':combined['external_validation_by_cell_group'],'scope':'Reporter-assay magnitude research model, with outside-test delivery protection. Raw projections/intervals are unvalidated diagnostics when the model withholds its authorized answer. No patient risk or causal gene output.'};write_json(OUT/'inference_registry.json',registry);return registry
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--input',required=True);parser.add_argument('--output',required=True);args=parser.parse_args();registry=json.loads((OUT/'inference_registry.json').read_text());rows=list(csv.DictReader(open(args.input)));assert rows,'No DNA inputs';valid=[];output=[];refs=[];alts=[];positions=[]
    for i,row in enumerate(rows):
        out={'id':row.get('id',str(i)),'cell_group':row.get('cell_group',''),'status':'unranked','unranked_reason':'','predicted_absolute_log2_effect':'','lower90':'','upper90':'','ranked_effect_magnitude':'','raw_model_projection':'','raw_interval_lower90':'','raw_interval_upper90':'','input_policy_eligible':0,'input_arm':registry['recommended_arm'],'AVI_available':''};output.append(out)
        if row.get('condition','').strip().lower() not in ['', 'baseline','untreated','unreported']:out['unranked_reason']='unsupported_stimulation_or_condition';continue
        a=row.get('reference','').upper();b=row.get('alternate','').upper();different=[j for j,(x,y) in enumerate(zip(a,b)) if x!=y]
        if len(a)!=len(b) or not 51<=len(a)<=1024 or set(a+b)-set('ACGT') or len(different)!=1:out['unranked_reason']='requires_one_actual_SNV_in_equal_unambiguous_51_to_1024bp_DNA_pairs';continue
        if row.get('cell_group') not in registry['cell_names']:out['unranked_reason']='unsupported_cell_group';continue
        valid.append(i);refs.append(a);alts.append(b);positions.append(different[0])
    if valid:
        unique={};inverse=[]
        for a,b,p in zip(refs,alts,positions):
            key=(a,b,p)
            if key not in unique:unique[key]=len(unique)
            inverse.append(unique[key])
        pairs=list(unique);encoded=Encoder().features([p[0] for p in pairs],[p[1] for p in pairs],[p[2] for p in pairs]);x=encoded[np.array(inverse,dtype='int32')];n=len(valid);d={'n':n,'x':x,'bank_index':np.arange(n,dtype='int32'),'unit':np.zeros(n,dtype='int8'),'cell_names':registry['cell_names'],'cell_line_names':registry['cell_line_names'],'family_names':registry['family_names'],'family_by_cell':np.array(registry['family_by_cell'],dtype='int16'),'n_context_flags':registry['n_context_flags'],'cell':np.array([registry['cell_names'].index(rows[i]['cell_group']) for i in valid],dtype='int16'),'line':np.array([registry['cell_line_names'].index(rows[i]['cell_line']) if rows[i].get('cell_line') in registry['cell_line_names'] else -1 for i in valid],dtype='int16')};avi=None
        if registry['recommended_arm']=='AVI_added':
            lookup={}
            with gzip.open(DATA/'AVI_inference_lookup.csv.gz','rt') as f:
                for r in csv.DictReader(f):lookup.setdefault(r['variant_GRCh38'],[]).append(r)
            avi=np.tile(np.array([np.nan,np.nan,0],dtype='float32'),(n,1))
            # When a genomic key is supplied, verify its substituted alleles
            # against the actual paired bases before using an AVI input.
            for j,i in enumerate(valid):
                key=rows[i].get('variant_GRCh38','');fields=key.split(':');p=positions[j];start=max(0,min(p-25,len(refs[j])-51))
                if len(fields)==4 and (fields[2],fields[3])==(refs[j][p],alts[j][p]):
                    matches=[r for r in lookup.get(key,[]) if r['reference51']==refs[j][start:start+51] and int(r['changed_position51'])==p-start]
                    if matches:avi[j]=[float(matches[0]['raw']),float(matches[0]['phred']),1]
        predictions=predict_ensemble(d,np.arange(n),registry['records'],avi);r=registry['records'][0]['record'];tx=Path(r['transform']);t=transform_load(tx if tx.is_absolute() else ROOT/tx);maxz=np.max(abs((x-t['mean'])/t['sd']),axis=1);cal=registry['calibration']
        for j,i in enumerate(valid):
            out=output[i];cell=str(d['cell'][j]);s=cal['cell_error_scale'].get(cell,cal['global_error_scale']);half=s*cal['worst_study_standardized_q90'];p=float(predictions[j]);lo=max(0.,p-half);hi=p+half;reasons=[]
            if cal['cell_training_study_support'].get(cell,0)<2:reasons.append('cell_group_supported_by_fewer_than_two_training_studies')
            if maxz[j]>12:reasons.append('sequence_features_outside_frozen_domain_limit')
            if hi-lo>cal['maximum_interval_width']:reasons.append('interval_too_wide')
            input_eligible=not reasons
            if not registry['cell_group_authorized_for_automated_research_ranking'].get(out['cell_group'],False):reasons.append('outside_tests_do_not_support_automatic_ranking_for_this_cell_group')
            answer=not reasons;out.update({'predicted_absolute_log2_effect':p if answer else '','lower90':lo if answer else '','upper90':hi if answer else '','raw_model_projection':p,'raw_interval_lower90':lo,'raw_interval_upper90':hi,'input_policy_eligible':int(input_eligible),'status':'research_candidate' if answer else 'unranked','unranked_reason':';'.join(reasons),'ranked_effect_magnitude':p if answer else '','AVI_available':int(avi[j,2]) if avi is not None else ''})
    with open(args.output,'w') as f:w=csv.DictWriter(f,fieldnames=list(output[0]));w.writeheader();w.writerows(output)
    print('RESEARCH MAGNITUDE SCORES WRITTEN',len(output),'ranked',sum(r['status']=='research_candidate' for r in output),flush=True)
if __name__=='__main__':main()

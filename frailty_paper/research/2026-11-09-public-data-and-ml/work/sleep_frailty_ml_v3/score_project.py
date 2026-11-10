"""Apply the locked final models to the33project substitutions in2contexts."""
import bisect,collections
from common import *
from raw_features import RawDNAInference
def main():
    registry_path=OUT/'final_model_registry.json';registry=json.loads(registry_path.read_text());lock=json.loads((OUT/'final_artifact_lock.json').read_text());assert lock['registry_sha256']==sha(registry_path);records=registry['records']
    for r in records.values():
        if r['artifact']:assert r['checkpoint_sha256']==sha(OUT/r['artifact'])
    input_path=OUT/'project_raw_inputs.json';receipt=json.loads((OUT/'project_input_receipt.json').read_text());assert receipt['output_sha256']==sha(input_path);inputs=json.loads(input_path.read_text());assert len(inputs)==66;engine=RawDNAInference(inputs);ix=np.arange(len(inputs));scores={n:engine.predict(records[n],ix,records) for n in set(registry['role_model_names'].values())};aa=engine.avi_features();training_coordinates=collections.defaultdict(list);training_variants=set()
    for stem in ['salomon','siraj']:
        for r in csv.DictReader((WORK/(stem+'_manifest_local.csv')).open()):
            if int(r['chromosome'])>=4:training_coordinates[int(r['chromosome'])].append(int(r['position_hg38']));training_variants.add(r['variant_GRCh38'])
    for c in training_coordinates:training_coordinates[c].sort()
    rows=[]
    for role,name in registry['role_model_names'].items():
        for i,r in enumerate(inputs):
            p=scores[name];c=int(r['chromosome']);pos=r['position_hg38'];arr=training_coordinates[c];j=bisect.bisect_left(arr,pos-1000);near=j<len(arr) and arr[j]<=pos+1000;conditional='conditional_scores' in p;variant_score=float(p['conditional_scores'][i,1] if conditional else p['overall'][i]);rows.append({'variant':r['rsid'],'variant_GRCh38':r['variant_GRCh38'],'assay_context':r['assay_context'],'model_role':role,'model_name':name,'overall_assay_call_score_uncalibrated':float(p['overall'][i]),'active_region_variant_score_uncalibrated':variant_score,'variant_score_definition':'Conditional active-region expert' if conditional else 'Overall call score; active evaluation uses a restricted subset','avi_raw':float(aa[i,0]) if aa[i,2] else '', 'avi_phred':float(aa[i,1]) if aa[i,2] else '', 'AVI_available':bool(aa[i,2]),'same_substitution_present_in_v3_training':r['variant_GRCh38'] in training_variants,'within_1000bp_of_v3_training_variant':near,'allele_direction':'Genome reference to specified alternate; not an original study risk/effect direction','interpretation':'Molecular experiment-prioritization hypothesis in this lab context. No clinical risk, causal gene or independent genetic-score validation.'})
    target=OUT/'project_final_model_predictions.csv'
    with target.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (OUT/'project_final_prediction_receipt.json').write_text(json.dumps({'variants':33,'assay_contexts':STUDIES,'rows':len(rows),'roles':registry['role_model_names'],'registry_sha256':sha(registry_path),'input_sha256':sha(input_path),'output_sha256':sha(target),'AVI_exact_allele_available_context_rows':int(aa[:,2].sum()),'AVI_source_receipts':engine.avi_source_receipts,'measured_accuracy_independently_established_for_candidates':False,'training_overlap_caution':'Rows explicitly flag exact-substitution and nearby training exposure. None are independent genetic-score validation.','scores_can_be_compared':'Within a fixed role/model and assay context; prevalence/detection differs between studies, so do not rank cell-type relevance from raw between-context values.','gene_assignment':'No causal gene assignment or directional disease-risk inference is made by this score.','source_provenance':'Original UCSC payload and exact target alleles verified in project_input_receipt.json'},indent=2)+'\n');print('Scored33project variants in2contexts using locked models;',len(rows),'role/context rows',flush=True)
if __name__=='__main__':main()

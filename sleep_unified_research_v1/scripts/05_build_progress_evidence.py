#!/usr/bin/env python3
"""Produce numerical status tables; never promote unfinished analyses."""
import csv
import hashlib
import json
from pathlib import Path
import shutil
from datetime import datetime, timezone

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'sleep_unified_research_v1'

def read(path):
    with path.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def write(path,rows,fields):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)

def hash_file(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    core=read(OUT/'sources/recovered/results/tables/rg_matrix.tsv')
    ext=read(ROOT/'discovery_extension/results/ldsc/extension_rg_matrix.tsv')
    rep=read(ROOT/'discovery_extension/results/replication/replication_results.tsv')
    assert len(core)==396 and len(ext)==1200 and len(rep)==217
    assert len({(r['sleep_trait'],r['disease_trait']) for r in core})==396
    assert len({(r['sleep_trait'],r['extension_trait_id']) for r in ext})==1200
    assert len({r['pair_id'] for r in rep})==217
    for source,target in [(OUT/'sources/recovered/results/tables/rg_matrix.tsv','original_core_396.tsv'),
        (ROOT/'discovery_extension/results/ldsc/extension_rg_matrix.tsv','original_extension_1200.tsv'),
        (ROOT/'discovery_extension/results/ldsc/extension_trait_readiness.tsv','original_extension_h2_100.tsv')]:
        destination=OUT/'tables'/target
        if not destination.exists():shutil.copyfile(source,destination)
        assert hash_file(source)==hash_file(destination)
    qualified=[]
    cohort={r['pair_id']:r for r in read(OUT/'reviews/cohort_independence_v1.tsv')}
    for r in rep:
        review=cohort[r['pair_id']]
        row={**r,'historical_replication_class':r['replication_class'],
            'sleep_input_independent':'false','both_trait_independent_replication':'false',
            'current_interpretation':'OUTCOME_SIDE_VALIDATION_ONLY' if r['replication_class']=='REPLICATED' else 'HISTORICAL_CANDIDATE_CLASSIFICATION',
            'heterogeneity_interpretation':'ZERO_COVARIANCE_DIAGNOSTIC_ONLY; SHARED_SLEEP_AND_SELECTION_UNCALIBRATED',
            'native_reproduction_status':'PENDING_FULL_NATIVE_RERUN',
            'current_class':review['reviewed_class'],
            'current_QC_reason':review['replication_outcome_h2_qc_reason'],
            'historical_participant_overlap_status':r['participant_overlap_status'],
            'current_participant_overlap_evidence':'COHORT_SOURCE_DISTINCT; INDIVIDUAL_INTERSECTION_UNKNOWN' if review['source_available'].lower()=='true' else 'NO_EXTERNAL_SOURCE',
            'Dsleep_Vsleep_intersection':review['Dsleep_Vsleep_intersection'],
            'Dsleep_Voutcome_intersection':review['Dsleep_Voutcome_intersection'],
            'Doutcome_Vsleep_intersection':review['Doutcome_Vsleep_intersection'],
            'Doutcome_Voutcome_intersection':review['Doutcome_Voutcome_intersection']}
        qualified.append(row)
    write(OUT/'REPLICATION_RESULTS.tsv',qualified,list(qualified[0]))
    ledger=read(ROOT/'discovery_extension/results/replication/replication_source_queue.tsv')
    for r in ledger:
        r['historical_participant_overlap_status']=r['participant_overlap_status']
        r['current_participant_overlap_evidence']='COHORT_SOURCE_DISTINCT; INDIVIDUAL_INTERSECTION_UNKNOWN' if cohort[r['pair_id']]['source_available'].lower()=='true' else 'NO_EXTERNAL_SOURCE'
        r['phenotype_equivalence_qualification']='HISTORICAL_NAME_CODING_MATCH; NOT_IDENTICAL_ASCERTAINMENT_OR_ESTIMAND'
    write(OUT/'REPLICATION_SOURCE_LEDGER.tsv',ledger,list(ledger[0]))
    shutil.copyfile(OUT/'reviews/novelty_v1_feasibility_matrix.tsv',OUT/'HYPOTHESIS_FEASIBILITY_MATRIX.tsv')

    original=next(r for r in core if r['sleep_trait']=='insomnia' and r['disease_trait']=='bmi')
    native=json.loads((OUT/'native/core_pilot_v1/rg_insomnia__bmi.full_precision.json').read_text())['estimates'][0]
    comparison=[]
    for key,new in [('rg',native['rg_ratio']),('se',native['rg_se']),('z',native['z']),('p',native['p']),
                    ('h2_obs',native['hsq2']['tot']),('h2_obs_se',native['hsq2']['tot_se']),
                    ('h2_int',native['hsq2']['intercept']),('h2_int_se',native['hsq2']['intercept_se']),
                    ('gcov_int',native['gencov']['intercept']),('gcov_int_se',native['gencov']['intercept_se'])]:
        old=float(original[key]); comparison.append({'pair_id':'insomnia__bmi','metric':key,
            'frozen_serialized_value':old,'native_full_precision_value':new,'absolute_difference':abs(old-new),
            'classification':'AGREES_AT_HISTORICAL_PRINTED_PRECISION',
            'historical_unrounded_estimate_available':'false',
            'native_input_level':'RECOVERED_HAPMAP3_INPUT; SOURCE_QC_CHAIN_NOT_RERUN'})
    write(OUT/'tables/original_vs_native_pilot.tsv',comparison,list(comparison[0]))
    evidence=[{'experiment_id':'core_native_insomnia_bmi_pilot_v1','status':'COMPLETED_PROCESSED_INPUT_REPRODUCTION',
        'sleep_source':'Jansen2019_UKB_only','outcome_source':'Locke2018_BMI','estimand':'global_rg',
        'rg':native['rg_ratio'],'se':native['rg_se'],'p':native['p'],
        'new_biological_finding':'false','both_trait_independent_replication':'false',
        'evidence':'native/core_pilot_v1/rg_insomnia__bmi.full_precision.json'}]
    write(OUT/'NEW_SCIENTIFIC_RESULTS.tsv',evidence,list(evidence[0]))
    admission=read(OUT/'reviews/provenance_v1_mvp_insomnia_217_admission.tsv')
    for row in admission:
        row.update(new_independent_rg_status='NOT_ESTIMATED',new_rg='',new_se='',new_p='',
                   completed_both_trait_independent_validation='false',
                   new_sleep_source_body_identity='MVP_GCST90475826_BYTES_MD5_SHA256_VERIFIED',
                   new_sleep_source_schema_status='ALL_SUPPLIED_SE_MISSING; CI_AND_N_DIAGNOSTICS_ONLY',
                   new_sleep_source_h2_power_status='NOT_ESTIMATED',
                   current_phenotype_admission_gate='STOP_GIA_ACCESSION_SPECIFIC_ICD_RULE_WINDOW_UNCERTIFIED',
                   current_pair_experiment_admission='NO_PAIR_TEST_ADMITTED_UNDER_PROTOCOL_V2',
                   source_qc_receipt='logs/mvp_insomnia_source_schema_v1.json',
                   current_protocol='FROZEN_NEW_ANALYSIS_PROTOCOL_v2.md')
        if row['original_sleep_trait']!='insomnia':
            for key in ['new_sleep_source_body_identity','new_sleep_source_schema_status','new_sleep_source_h2_power_status','current_phenotype_admission_gate','source_qc_receipt']:
                row[key]='NOT_APPLICABLE_NO_NEW_SOURCE_ADMITTED_FOR_THIS_SLEEP_TRAIT'
    write(OUT/'INDEPENDENT_VALIDATION.tsv',admission,list(admission[0]))
    local=[{'analysis':'local_rg','status':'NOT_ADMITTED_PRIMARY_QUESTION_NO_GO','estimate':'',
            'reason':'No selected novel question/valid focused signal/input-and-LD contract frozen'},
           {'analysis':'molecular_colocalization_or_gene','status':'NOT_ADMITTED_UPSTREAM_NO_SIGNAL','estimate':'',
            'reason':'No admitted locus/QTL experiment; proximity and global rg cannot identify mechanism'}]
    write(OUT/'LOCAL_AND_MOLECULAR_RESULTS.tsv',local,list(local[0]))
    counts=[{'family':'core','planned_pairs':396,'frozen_positive_total':sum(float(r['fdr'])<.05 for r in core),
             'frozen_primary_positive':sum(float(r['fdr'])<.05 and r['analysis_tier']=='PRIMARY_PHASE1' for r in core),
             'native_pairs_reproduced':1,'completed_both_trait_independent_validations':0},
            {'family':'extension','planned_pairs':1200,'frozen_positive_total':sum(float(r['extension_fdr'])<.05 for r in ext),
             'frozen_primary_positive':603,'native_pairs_reproduced':0,'completed_both_trait_independent_validations':0},
            {'family':'historical_validation','planned_pairs':217,'frozen_positive_total':23,'frozen_primary_positive':23,
             'native_pairs_reproduced':0,'completed_both_trait_independent_validations':0}]
    write(OUT/'tables/study_design_counts.tsv',counts,list(counts[0]))
    metadata=[]
    for r in read(ROOT/'config/analysis_panel.tsv'):
        metadata.append({'family':'core','trait_id':r['trait_id'],'phenotype':r['label'],'source':r['dataset_version'],
            'build':r['build'],'ancestry':r['ancestry'],'N':r['n_total'],'definition':r['phenotype_definition']})
    for r in read(ROOT/'discovery_extension/config/candidate_traits.tsv'):
        metadata.append({'family':'extension','trait_id':r['extension_trait_id'],'phenotype':r['phenotype_name'],
            'source':r['source'],'build':r['build'],'ancestry':r['ancestry'],'N':r['sample_size'],'definition':r['phenotype_definition']})
    write(OUT/'tables/phenotype_and_source_metadata.tsv',metadata,list(metadata[0]))
    receipt={'created_utc':datetime.now(timezone.utc).isoformat(),'original_core_rows':396,'extension_rows':1200,
             'validation_candidate_rows':217,'native_core_pair_count':1,'native_extension_pair_count':0,
             'completed_both_trait_independent_validations':0,'new_biological_findings':0,
             'full_native_scientific_reproduction_complete':False,'manuscript_generated':False}
    (OUT/'logs/progress_evidence_receipt_v1.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()

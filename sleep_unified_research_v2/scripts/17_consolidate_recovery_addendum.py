#!/usr/bin/env python3
"""Consolidate completed recovery evidence without changing sealed v1."""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1]
V1 = PACKAGE.parent/'sleep_unified_research_v1'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))

def load(relative):
    return json.loads((PACKAGE/relative).read_text())

def save(path,value):
    with path.open('x') as f:
        json.dump(value,f,indent=2)
        f.write('\n')

def table(path,rows):
    with path.open('x',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def main():
    receipt=PACKAGE/'logs/recovery_consolidation_receipt_v2.json'
    if receipt.exists():
        raise SystemExit('PRIOR_CONSOLIDATION_PRESERVED')
    old=read(V1/'tables/native_input_hash_checks.tsv')
    panel={r['trait_id']:r for r in read(PACKAGE.parent/'config/analysis_panel.tsv')}
    schemas={r['trait_id']:r for r in read(PACKAGE.parent/'config/gwas_schemas.tsv')}
    sources={trait:r for r in read(PACKAGE.parent/'config/public_gwas_sources.tsv') for trait in r['trait_ids'].split(',')}
    def single_trait_metadata(trait):
        p=panel[trait];s=schemas[trait];r=sources[trait]
        return dict(source_release=p['dataset_version'],source_id=p['source_id'],
                    source_page_url=r['source_page_url'],genome_build=p['build'],
                    source_build_verification_status=r['build_status'],ancestry=p['ancestry'],
                    phenotype_identity=p['phenotype_definition'],effect_coding=s['effect_convention'],
                    permitted_use='RESEARCH_ONLY;SOURCE_TERMS_AND_ATTRIBUTION_REQUIRED;NO_RAW_REDISTRIBUTION',
                    required_downstream_tasks='COMPLETE_HARMONIZATION_MUNGING_AND_NATIVE_ESTIMATION_REPLAY')
    def metadata(trait_ids):
        traits=trait_ids.split(',')
        rows=[single_trait_metadata(trait) for trait in traits]
        return {key:rows[0][key] if len({r[key] for r in rows})==1 else ' | '.join(
                    trait+': '+r[key] for trait,r in zip(traits,rows)) for key in rows[0]}
    missing=[r for r in old if r['status']=='MISSING']
    recovered={}
    for r in load('logs/additional_ssd_sources_verification_v2.json')['rows']:
        if r['status']=='EXACT_EXPECTED_SHA256_READ_ONLY_SOURCE_RECOVERED':
            recovered.setdefault(r['actual_sha256'],dict(path=r['path'],bytes=r['bytes'],actual_sha256=r['actual_sha256'],method='READ_ONLY_EXISTING_SSD_FILE'))
    for r in load('logs/direct_tsv_source_recovery_receipt_v2.json')['results']:
        if r['status']=='EXACT_SOURCE_BYTES_RESTORED_FROM_HISTORICAL_GZIP':
            recovered[r['actual_source_sha256']]=dict(path=r['path'],bytes=r['actual_source_bytes'],actual_sha256=r['actual_source_sha256'],method='EXACT_DECOMPRESSED_SOURCE_COPY_IN_NEW_SSD_NAMESPACE')
    acquired=load('logs/exact_lipid_ranged_acquisition_receipt_v2.json')
    if not acquired['all_exact_expected_hashes'] or len(acquired['sources'])!=3:
        raise SystemExit('ALL_THREE_EXACT_LIPID_SOURCES_REQUIRED')
    for r in acquired['sources']:
        if r['status']!='EXACT_HISTORICAL_SOURCE_REACQUIRED' or r['actual_sha256']!=r['expected_sha256']:
            raise SystemExit('EXACT_SOURCE_GATE_FAILED')
        recovered[r['actual_sha256']]=dict(path=r['path'],bytes=r['actual_bytes'],actual_sha256=r['actual_sha256'],method='EXACT_HISTORICAL_PUBLIC_FILE_REACQUIRED')
    updates=[]
    for r in missing:
        match=recovered.get(r['expected_sha256'])
        updates.append(dict(kind=r['kind'],trait_id=r['trait_id'],original_missing_path=r['path'],
                            expected_sha256=r['expected_sha256'],resolved_path=match['path'] if match else '',
                            resolved_bytes=match['bytes'] if match else '',actual_sha256=match['actual_sha256'] if match else '',
                            recovery_method=match['method'] if match else 'UNRECOVERED_ORIGINAL_SOURCE_CONTAINER',
                            current_status='EXACT_EXPECTED_SOURCE_BYTES_AVAILABLE_SEPARATE_LOCATION' if match else 'UNRECOVERED',
                            original_missing_path_replaced=False,source_preprocessing_or_estimation_implied=False,
                            **metadata(r['trait_id'])))
    table(PACKAGE/'tables/SOURCE_RECOVERY_UPDATES.tsv',updates)
    current=[]
    for r in old:
        if r['kind'] not in {'core_raw','core_source_archive'}:
            continue
        match=recovered.get(r['expected_sha256']) if r['status']=='MISSING' else None
        current.append(dict(kind=r['kind'],trait_id=r['trait_id'],expected_sha256=r['expected_sha256'],
                            current_path=match['path'] if match else r['path'],
                            current_bytes=match['bytes'] if match else r['bytes'],
                            current_sha256=match['actual_sha256'] if match else r['actual_sha256'],
                            current_status='EXACT_HISTORICAL_HASH_AVAILABLE' if match or r['status']=='MATCH_EXPECTED_SHA256' else 'MISSING',
                            proof='tables/SOURCE_RECOVERY_UPDATES.tsv' if match else 'sealed_v1/tables/native_input_hash_checks.tsv',
                            **metadata(r['trait_id'])))
    table(PACKAGE/'tables/CURRENT_CORE_SOURCE_STATUS.tsv',current)
    pref=[]
    for trait in ['hdl','ldl','triglycerides']:
        r=load('logs/'+trait+'_native_prefilter_receipt_v2.json')
        if r['status']!='NATIVE_PREFILTER_EXACT_HISTORICAL_BYTES_AND_COUNTS':
            raise SystemExit('ALL_THREE_NATIVE_PREFILTER_COMPARISONS_REQUIRED')
        c=r['comparison']
        pref.append(dict(trait_id=trait,historical_source_rows=c['source_rows'][0],native_source_rows=c['source_rows'][1],
                         historical_retained_rows=c['retained_rows'][0],native_retained_rows=c['retained_rows'][1],
                         historical_output_sha256=c['gzip_SHA256'][0],native_output_sha256=c['gzip_SHA256'][1],
                         status=r['status'],elapsed_seconds=r['elapsed_seconds'],worker_peak_observed_RSS_bytes=r['worker_peak_observed_RSS_bytes'],
                         before_after_input_hashes_unchanged=r['input_hashes_unchanged'],complete_GWAS_chain_reproduced=False,
                         receipt='logs/'+trait+'_native_prefilter_receipt_v2.json'))
    table(PACKAGE/'tables/NATIVE_LIPID_PREFILTER_COMPARISON.tsv',pref)
    phases=[
        (0,'QUALIFIED_RECOVERY_COMPLETE','45/45 raw hashes;40/43 source containers;3 originals unavailable','tables/CURRENT_CORE_SOURCE_STATUS.tsv'),
        (1,'PARTIAL','18/18 checkpoint artifacts;1 native core rg control;3 native prefilters;full45-trait QC/replay pending','sealed_v1/CORE_NATIVE_REPRODUCTION.md;tables/NATIVE_LIPID_PREFILTER_COMPARISON.tsv'),
        (2,'PARTIAL','100h2 and1200rg printed-log agreement;603BH;full preprocessing/native fits pending','sealed_v1/EXTENSION_NATIVE_REPRODUCTION.md'),
        (3,'PARTIAL','217classified;23positive18concordant17QC159unavailable;0both-trait-independent successes','sealed_v1/REPLICATION_RESULTS.tsv;FROZEN_NEW_ANALYSIS_PROTOCOL_v3.md'),
        (4,'METHOD_BLOCKED','zero-covariance diagnostics qualified;aligned genomic blocks and calibrated ratio covariance missing','sealed_v1/SHARED_ESTIMATOR_COVARIANCE.md'),
        (5,'DECISION_COMPLETE_NO_GO','No eligible novel primary question;future lead requires data/power/covariance/independent test','sealed_v1/HYPOTHESIS_FEASIBILITY_MATRIX.tsv'),
        (6,'NO_EXPERIMENT_ADMITTED','Ten conditional clinical transport candidates remain gated;0new pair outcomes','FROZEN_NEW_ANALYSIS_PROTOCOL_v3.md'),
        (7,'NOT_ADMITTED','No justified eligible local experiment','sealed_v1/LOCAL_AND_MOLECULAR_RESULTS.tsv'),
        (8,'NOT_ADMITTED','No genomic signal admitted for molecular/clinical follow-up','sealed_v1/LOCAL_AND_MOLECULAR_RESULTS.tsv'),
        (9,'PARTIAL','Verified frozen tables and7traceable figures;complete native comparisons unavailable','sealed_v1/FINAL_EVIDENCE_HANDOFF.md'),
        (10,'AUTOMATED_REVIEW_COMPLETE_QUALIFIED','Seven separate v1 roles;additional source and recovery/numerical reviews;human questions unresolved','sealed_v1/ADVERSARIAL_REVIEW.md;source_definition_review/;reviews/'),
        (11,'DECISION_COMPLETE_NOT_READY','No manuscript-authoring clearance;scientific completeness and SLEEP novelty requirements unmet','SCIENTIFIC_COMPLETION_V2.md'),
        (12,'ADDENDUM_PREPARED_FOR_VALIDATION','V1 sealed;v2 technical reports/receipts/tests;commit/export proof supplied separately','README.md')]
    table(PACKAGE/'tables/PHASE_COMPLETION_V2.tsv',[dict(phase=p,status=s,evidence_summary=e,proof=f) for p,s,e,f in phases])
    status={'recorded_utc':datetime.now(timezone.utc).isoformat(),
            'initial_missing_ledger_entries':len(missing),'entries_resolved':sum(r['current_status']!='UNRECOVERED' for r in updates),
            'initial_unique_missing_hashes':len({r['expected_sha256'] for r in missing}),
            'unique_source_hashes_recovered':len(recovered),'remaining_original_containers':[r['trait_id'] for r in updates if r['current_status']=='UNRECOVERED'],
            'raw_sources_exact_available':sum(r['kind']=='core_raw' and r['current_status']=='EXACT_HISTORICAL_HASH_AVAILABLE' for r in current),
            'source_containers_exact_available':sum(r['kind']=='core_source_archive' and r['current_status']=='EXACT_HISTORICAL_HASH_AVAILABLE' for r in current),
            'native_prefilters_exact':len(pref),'native_core_rg_controls_total_including_v1':1,
            'complete_core_396_native_reproduction':False,'complete_extension_1200_native_reproduction':False,
            'complete_validation_41_native_reproduction':False,'new_clinical_pair_outcomes_examined':0,
            'both_trait_independent_successes':0,'clinical_definition_public_gate_resolved':True,
            'new_primary_question_admitted':False,'corrected_shared_covariance_calibrated':False,
            'local_or_molecular_experiments_admitted':False,'HUMAN_MANUSCRIPT_AUTHORING_READY':False,
            'full_native_resource_preflight':load('logs/native_resource_preflight_v2.json'),
            'original_SSD_and_sealed_v1_outputs_modified':False}
    save(PACKAGE/'tables/CURRENT_EVIDENCE_STATUS_V2.json',status)
    table_paths=['SOURCE_RECOVERY_UPDATES.tsv','CURRENT_CORE_SOURCE_STATUS.tsv','NATIVE_LIPID_PREFILTER_COMPARISON.tsv','PHASE_COMPLETION_V2.tsv','CURRENT_EVIDENCE_STATUS_V2.json']
    save(receipt,{'completed_utc':datetime.now(timezone.utc).isoformat(),'script_sha256':sha(Path(__file__)),
                  'tables_sha256':{p:sha(PACKAGE/'tables'/p) for p in table_paths},'current_status':status})
    print(json.dumps(status,indent=2))

if __name__=='__main__':
    main()

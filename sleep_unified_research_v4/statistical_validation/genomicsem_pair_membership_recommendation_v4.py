#!/usr/bin/env python3
"""Result-free membership recommendation derived from the frozen native plan.

Never uses rg/SE/P values or reads result tables, never estimates covariance, and retains all
217 candidate rows in the original manifest. Writes distinct v4 metadata only.
"""
from pathlib import Path
import csv
import json
import hashlib
from datetime import datetime, timezone

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'sleep_unified_research_v4/statistical_validation'


def main():
    plan_path=ROOT/'sleep_unified_research_v4/manifests/ssd_native_execution_plan_v4_3.json'
    manifest_path=ROOT/'discovery_extension/config/replication_manifest.tsv'
    plan=json.loads(plan_path.read_text())
    manifest=list(csv.DictReader(manifest_path.open(),delimiter='\t'))
    assert len(manifest)==217
    inputs={Path(r['path']).name:r for r in plan['inputs_verified']}
    assert len(inputs)==158
    rows=[]
    for job in plan['jobs']:
        if job['stage']!='validation' or job['kind']!='rg':continue
        sleep_file=Path(job['inputs'][0]);sleep=sleep_file.name.removesuffix('.sumstats.gz')
        for path in job['inputs'][1:]:
            val_file=Path(path);source=val_file.name.removesuffix('.sumstats.gz')
            matched=[r for r in manifest if r['sleep_trait']==sleep and r['replication_source_id']==source]
            assert len(matched)==1
            r=matched[0]
            disc=inputs[r['extension_trait_id']+'.sumstats.gz'];slp=inputs[sleep_file.name];val=inputs[val_file.name]
            assert all(x['match'] for x in [disc,slp,val])
            rows.append({'member_order':len(rows)+1,'native_validation_job_id':job['job_id'],
                'historical_pair_id':r['pair_id'],'sleep_trait':sleep,'discovery_outcome_id':r['extension_trait_id'],
                'validation_source_id':source,'phenotype_name':r['external_phenotype_name'],
                'sleep_input':slp['path'],'sleep_sha256':slp['actual_sha256'],
                'discovery_outcome_input':disc['path'],'discovery_outcome_sha256':disc['actual_sha256'],
                'validation_outcome_input':val['path'],'validation_outcome_sha256':val['actual_sha256'],
                'membership_rule':'EXACT_FROZEN_VALIDATION_RG_PLAN; no rg/SE/P selection',
                'new_covariance_status':'NOT_ESTIMATED; all common-boundary/native-calibration gates pending',
                'interpretation_if_eventually_estimated':'RETROSPECTIVE_SHARED_SLEEP_SOURCE_COMPARISON; original217 family retained'})
    assert len(rows)==41 and len({r['historical_pair_id'] for r in rows})==41
    target=OUT/'genomicsem_covariance_pair_members_v4.tsv'
    with target.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    controls=[]
    # Existing native control identity is the selection reason, not its rg or P.
    for name,sleep,outcome,reason in [
        ('CORE_INSOMNIA_BMI_DUPLICATE','insomnia','bmi','Existing separately verified native control; exact duplicate must give zero contrast'),
        ('FIRST_FROZEN_VALIDATION_PAIR_DUPLICATE',rows[0]['sleep_trait'],rows[0]['validation_source_id'],'First member of unchanged native validation plan; no new outcome-based selection')]:
        slp=inputs[sleep+'.sumstats.gz'];other=inputs[outcome+'.sumstats.gz']
        controls.append({'control_id':name,'sleep_trait':sleep,'outcome_id':outcome,
            'input1':slp['path'],'input1_sha256':slp['actual_sha256'],
            'input2':other['path'],'input2_sha256':other['actual_sha256'],
            'selection_reason':reason,'arms':'Pinned pair one-step and explicitly configured two-step in separately frozen calibration; compare duplicates within each arm',
            'required_checks':'Duplicate points/deletes/marginals; covariance equals marginal variance; paired delta and SE zero; default-vs-common200 point consistency within the same arm',
            'status':'RECOMMENDED_NOT_EXECUTED; prototype native guard remains active'})
    control_path=OUT/'genomicsem_native_calibration_controls_v4.tsv'
    with control_path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(controls[0]),delimiter='\t');w.writeheader();w.writerows(controls)
    receipt={'utc':datetime.now(timezone.utc).isoformat(),'original_candidate_count':217,
        'historical_native_comparison_members':41,'controls':len(controls),
        'new_result_tables_read':False,'new_estimates_computed':False,'new_primary_family_admitted':False,
        'original_manifest_modified':False,
        'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [plan_path,manifest_path]},
        'output_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [target,control_path]},
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (OUT/'genomicsem_pair_membership_recommendation_receipt_v4.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'candidate_count':217,'member_count':41,'controls':2,'new_covariance':'NOT_ESTIMATED'},indent=2))


if __name__=='__main__':main()

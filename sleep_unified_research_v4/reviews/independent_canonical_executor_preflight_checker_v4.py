#!/usr/bin/env python3
"""Read-only plan/dependency/control/prerequisite verification; no imports/fits."""
from collections import Counter
import datetime
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
P=ROOT/'sleep_unified_research_v4'
SEEN={}


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):
            h.update(b)
    return h.hexdigest()


def read(path):
    SEEN[str(path)]={'sha256':sha(path),'bytes':path.stat().st_size}
    return path.read_text()


def main():
    read(Path(__file__))
    pp=P/'manifests/native_canonical_calibration_plan_v4.json'
    plan=json.loads(read(pp));assert sha(pp)=='18baccdd5cea6ccdeceee6ec3ad386dcfb183693be89b908e730ccd4d0eb5d62'
    admission=json.loads(read(P/'manifests/native_canonical_calibration_admission_v4.json'))
    assert admission['execution_admitted'] is False and admission['plan_sha256']==sha(pp)
    assert admission['allow_41_covariance_outcomes'] is False and admission['allow_calibrated_biological_p_values'] is False
    assert plan['output_limit_bytes']==8*1024**3 and plan['worker_count']==plan['blas_threads']==1
    assert plan['internal_floor_bytes']==3*1024**3 and plan['ssd_floor_bytes']==5*1024**3
    assert plan['worker_rss_limit_bytes']==2*1024**3 and plan['per_worker_seconds_limit']==7200
    assert len(plan['jobs'])==16 and len({j['job_id'] for j in plan['jobs']})==16
    observed={(j['control_id'],j['method'],j['partition'],j['duplicate']) for j in plan['jobs']}
    expected={(c,m,p,d) for c in plan['allowed_control_ids'] for m in ('one_step','two_step30') for p in ('stock_default200','canonical200') for d in ('A','B')}
    assert observed==expected
    for j in plan['jobs']:
        assert len(j['inputs'])==2 and j['input_sha256']=={s:plan['inputs_current_sha256'][s] for s in j['inputs']}
        assert Path(j['output_dir']).parent==Path(plan['ssd_output_root'])
        opts=j['ldsc_args'];assert opts[opts.index('--rg')+1]==','.join(j['inputs'])
        assert opts[opts.index('--n-blocks')+1]=='200' and '--print-delete-vals' in opts
        assert ('--two-step' in opts)==(j['method']=='two_step30')
        if j['method']=='two_step30':assert opts[opts.index('--two-step')+1]=='30'
    dependency_actual={}
    for s,h in plan['dependencies_sha256'].items():
        q=Path(s);SEEN[s]={'sha256':sha(q),'bytes':q.stat().st_size};dependency_actual[s]=SEEN[s]['sha256'];assert dependency_actual[s]==h
    assert len(dependency_actual)==71
    historical=json.loads(read(Path(plan['historical_plan_path'])))
    assert len(historical['jobs'])==190 and sha(Path(plan['historical_plan_path']))==plan['historical_plan_sha256']
    moved=dict(historical['dependencies_sha256'])
    old=str(ROOT/'sleep_unified_research_v1/scripts/native_ldsc_capture.py')
    new=str(Path(historical['ssd_support_package'])/'scripts/native_ldsc_capture.py')
    moved[new]=moved.pop(old)
    assert sha(Path(old))==sha(Path(new))==moved[new]
    counts=Counter();literal_failures=[];mapped_passes=[];missing=[]
    for j in historical['jobs']:
        q=Path(historical['ssd_support_package'])/'native'/(j['stage']+'_reproduction_v1')/(j['job_id']+'.execution_receipt.json')
        if not q.exists():missing.append(j['job_id']);continue
        r=json.loads(read(q));counts[j['stage']]+=1
        if r['dependency_sha256_before']!=historical['dependencies_sha256'] or r['dependency_sha256_after']!=historical['dependencies_sha256']:
            literal_failures.append(j['job_id'])
        if r['dependency_sha256_before']==moved and r['dependency_sha256_after']==moved:
            mapped_passes.append(j['job_id'])
    snapshots={}
    for stem in ['42_prepare_native_canonical_calibration.py','canonical_calibration_common.py','canonical_native_adapter.py','canonical_calibration_capture.py','canonical_calibration_arithmetic.py','canonical_calibration_fixture.py']:
        q=P/'scripts'/stem;snapshots[str(q)]={'sha256':sha(q),'content':read(q)}
    snapshot=P/'reviews/independent_canonical_executor_original_sources_v4.json';assert not snapshot.exists()
    snapshot.write_text(json.dumps(snapshots,indent=2)+'\n')
    unchanged=all(sha(Path(s))==m['sha256'] for s,m in SEEN.items());assert unchanged
    result={'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'verdict':'PREPARATION_ONLY_PASS_EXECUTION_BINDING_AND_CLEANUP_AMENDMENT_REQUIRED',
            'plan_sha256':sha(pp),'admission_false_verified':True,'controls16_matrix_verified':True,'dependency71_hashes_match':True,
            '8GiB_operational_amendment_explicit':True,'historical_receipts_present_by_stage':dict(counts),'historical_receipts_missing':len(missing),
            'literal_original_dependency_map_gate_failures':len(literal_failures),'explicit_relocation_map_gate_passes':len(mapped_passes),
            'relocated_key':{'original':old,'actual':new,'sha256':moved[new]},
            'all190_completion_or_numerical_pass_claimed':False,'new_native_controls_executed':False,
            'full_point_mask_preservation':'Adapter retains original _rg rows,args,filter and regression constructors; only IRWLS separators change in canonical arm. Cross-arm order/mask/weighted-row proof should be explicit. Weighted arrays may differ in last bits due to block-sum reduction and two-step intercept propagation.',
            'one_step_default_alarm_withdrawn':{'reason':'estimate_rg deepcopies and normalizes intercept_h2 with _split_or_none before its default condition. None becomes [None,None], so subsequent intercept_h2 is None is false. One-step remains one-step. No default-suppression change is needed.','reviewer_error_preserved':True,'recommended_noninvasive_gate':'Assert effective method and free intercept before regression; fixture should exercise estimate_rg normalization.'},
            'findings':['Historical receipt dependency key must map exactly the original wrapper path to relocated support wrapper path.','Unhandled teardown can skip failure receipt and shared context explicitly unlocks despite survivors.','Add final post-teardown resource/deadline/output/admission gates and reject unexpected descendant groups even after parent return0.','Historical completion must verify output/full-capture hashes and stage monitor plan SHA in addition to flags/maps.','Cross-arm final SNP order/masks/Nbar/M exact checks and weighted-row numerical comparisons with declared tolerances should substantiate full-point/mask preservation.'],
            '8GiB_resource_limit_scientific_scope':'Explicit separate operational cap only; no biological threshold or family amendment. Parent root review still required before admission.',
            'realistic_LD_sampling_calibration_pass':False,'signed_genotype_LD_asset_verified':False,'allow41_covariance_outcomes':False,'new_calibrated_biological_pvalues_admitted':False,
            'source_snapshot_sha256':sha(snapshot),'all_consumed_hashes_unchanged':unchanged,'inputs':SEEN,'candidate_code_imported_or_executed':False,'GWAS_outcome_columns_opened':False,'hash_buffer_bytes':65536}
    out=P/'reviews/independent_canonical_executor_preflight_receipt_v4.json';assert not out.exists();out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['inputs','findings']},indent=2))


if __name__=='__main__':main()

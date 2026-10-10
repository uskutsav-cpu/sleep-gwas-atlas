#!/usr/bin/env python3
"""Scalar/delete arithmetic for the single qualified h2; no LDSC import/fit."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def scalar_deletes(value):
    result=[]
    for x in value:
        if isinstance(x,list):
            if len(x)!=1:raise RuntimeError('EXPECTED_ONE_COLUMN_DELETE_ARRAY')
            x=x[0]
        if not isinstance(x,(int,float)) or not math.isfinite(x):raise RuntimeError('NONFINITE_DELETE_STATISTIC')
        result.append(x)
    if len(result)!=200:raise RuntimeError('EXPECTED_200_DELETIONS')
    return result


def centered_se(value):
    mean=math.fsum(value)/len(value)
    return math.sqrt((len(value)-1)/len(value)*math.fsum((x-mean)**2 for x in value))


def verify(plan_path,expected_hash,out):
    if sha(plan_path)!=expected_hash:raise RuntimeError('FROZEN_H2_PLAN_CHANGED')
    plan=json.loads(plan_path.read_text())
    if plan['stage']!='observed_h2' or plan['scope']!='FINNGEN_INSOMNIA_SOURCE_FEASIBILITY_ONLY' or plan['new_rg_commands']!=0:
        raise RuntimeError('SINGLE_QUALIFIED_OBSERVED_H2_SCOPE_REQUIRED')
    stage=Path(plan['stage_receipt']);completed=json.loads(stage.read_text())
    if completed['status']!='QUALIFIED_FEASIBILITY_STAGE_COMPLETE_VERIFIED' or completed['plan_sha256']!=expected_hash or not completed['owned_cleanup_verified']:
        raise RuntimeError('H2_STAGE_NOT_COMPLETE_AND_REAPED')
    for p,h in completed['result_sha256'].items():
        if sha(p)!=h:raise RuntimeError('NATIVE_PROOF_CHANGED')
    job=plan['jobs'][0];capture=Path(job['result_receipt']);r=json.loads(capture.read_text())
    if len(r['estimates'])!=1 or r['plan_sha256']!=expected_hash or r['job_id']!=job['job_id']:
        raise RuntimeError('H2_CAPTURE_IDENTITY_OR_CARDINALITY_CHANGED')
    if r['input_sha256_before']!=plan['input_sha256'] or r['input_sha256_after']!=plan['input_sha256']:
        raise RuntimeError('SOURCE_DERIVATIVE_CHANGED')
    x=r['estimates'][0]
    if x['input']!=plan['derivative'] or x['n_blocks']!=200 or x['constrain_intercept']:
        raise RuntimeError('SINGLE_INPUT_OR_FREE_INTERCEPT_CHANGED')
    for k in ['tot','tot_se','intercept','intercept_se','mean_chisq','lambda_gc']:
        if not isinstance(x[k],(int,float)) or not math.isfinite(x[k]):raise RuntimeError('NONFINITE_H2_CAPTURE')
    if x['tot_se']<=0 or x['intercept_se']<=0:raise RuntimeError('NONPOSITIVE_STANDARD_ERROR')
    checks={}
    for deletion,field in [('tot_delete_values','tot_se'),('intercept_delete_values','intercept_se')]:
        value=centered_se(scalar_deletes(x[deletion]))
        checks[field]=dict(independent=value,native=x[field],fixed_tolerance_pass=math.isclose(value,x[field],rel_tol=1e-12,abs_tol=1e-15))
    z=x['tot']/x['tot_se']
    diagnostic=dict(trait='finngen_R13_F5_INSOMNIA',source='R13_generation_1777989563097164',
                    scale='OBSERVED_ONLY_WITH_ASSUMED_CONSTANT_EFFECTIVE_N',h2=x['tot'],h2_SE=x['tot_se'],h2_Z=z,
                    h2_CI_lower=x['tot']-1.959963984540054*x['tot_se'],h2_CI_upper=x['tot']+1.959963984540054*x['tot_se'],
                    intercept=x['intercept'],intercept_SE=x['intercept_se'],lambda_GC=x['lambda_gc'],mean_chi_square=x['mean_chisq'],
                    native_SNP_count=x['n_snp'],assumed_effective_N=plan['assumed_effective_N'],
                    necessary_h2_Z_ge4=z>=4,necessary_intercept_le1_2=x['intercept']<=1.2,
                    source_INFO_gate_verified=False,per_variant_N_verified=False,scientific_source_admitted=False,
                    fully_independent_two_trait_replication=False,new_rg_estimates=0,
                    diagnostic_scope='Source feasibility only; threshold conditions do not resolve INFO/N, ascertainment, overlap, clinical construct or power.')
    result=dict(plan_sha256=expected_hash,stage_receipt_sha256=sha(stage),full_precision_capture_sha256=sha(capture),
                verifier_sha256=sha(Path(__file__)),fixed_tolerances=dict(relative=1e-12,absolute=1e-15),
                arithmetic_checks=checks,arithmetic_pass=all(v['fixed_tolerance_pass'] for v in checks.values()),diagnostic=diagnostic,
                warning_lines=r['warnings_and_errors'],all_200_deletions_preserved_in_capture=True,
                common_cross_trait_boundaries_established=False,original_families_changed=False)
    with Path(out).open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    table=Path(out).with_suffix('.tsv')
    with table.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(diagnostic),delimiter='\t',lineterminator='\n');w.writeheader();w.writerow(diagnostic)
    if not result['arithmetic_pass']:raise RuntimeError('H2_DELETE_ARITHMETIC_FAILURE_PRESERVED')
    print(json.dumps(dict(arithmetic_pass=True,source_admitted=False,new_rg_estimates=0)))


def main():
    p=argparse.ArgumentParser();p.add_argument('--plan',type=Path,required=True);p.add_argument('--plan-sha256',required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();verify(a.plan,a.plan_sha256,a.out)


if __name__=='__main__':main()

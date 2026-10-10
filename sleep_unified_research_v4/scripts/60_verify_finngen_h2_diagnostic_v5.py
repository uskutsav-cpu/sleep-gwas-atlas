#!/usr/bin/env python3
"""Independent scalar arithmetic and stock-delete closure for one diagnostic h2."""
import argparse
import csv
import hashlib
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import re
from terminal_commit_common_v2 import TerminalCommit, require_committed

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def regular(path):
    p=Path(path)
    if p.is_symlink() or not p.is_file() or any(parent.is_symlink() for parent in p.parents):
        raise RuntimeError('EXACT_REGULAR_DIAGNOSTIC_EVIDENCE_REQUIRED: '+str(p))
    return p

def regular_hashes(mapping):
    if not isinstance(mapping,dict) or not mapping:raise RuntimeError('NONEMPTY_FROZEN_RESULT_MAP_REQUIRED')
    for path,digest in mapping.items():
        if not isinstance(digest,str) or not re.fullmatch('[0-9a-f]{64}',digest) or sha(regular(path))!=digest:
            raise RuntimeError('FROZEN_DIAGNOSTIC_EVIDENCE_CHANGED: '+path)

def scalar_deletes(value):
    result=[]
    for x in value:
        if isinstance(x,list):
            if len(x)!=1:raise RuntimeError('EXPECTED_ONE_COLUMN_DELETE_ARRAY')
            x=x[0]
        if type(x) not in (int,float) or not math.isfinite(x):raise RuntimeError('NONFINITE_DELETE_STATISTIC')
        result.append(x)
    if len(result)!=200:raise RuntimeError('EXPECTED_200_DELETIONS')
    return result

def stock_scalar_deletes(path):
    with regular(path).open() as f:
        lines=[line.split() for line in f if line.strip()]
    if len(lines)!=200 or any(len(row)!=1 for row in lines):
        raise RuntimeError('STOCK_200_SINGLE_COLUMN_DELETIONS_REQUIRED')
    return scalar_deletes([float(row[0]) for row in lines])

def centered_se(value):
    mean=math.fsum(value)/len(value)
    return math.sqrt((len(value)-1)/len(value)*math.fsum((x-mean)**2 for x in value))

def report_terminal_paths(out):
    return Path(str(out)+'.pending.json'),Path(str(out)+'.terminal.json')

def require_verified_report(plan_path,expected_hash,out):
    plan_path=regular(plan_path)
    if sha(plan_path)!=expected_hash:raise RuntimeError('FROZEN_REPORT_INPUT_PLAN_CHANGED')
    plan=json.loads(plan_path.read_text());out=regular(out);table=regular(out.with_suffix('.tsv'))
    r=json.loads(out.read_text());stage=regular(plan['stage_receipt']);stage_hash=sha(stage)
    native=json.loads(stage.read_text());executor=regular(plan['executor_path'])
    if len(plan['jobs'])!=1 or plan['stage']!='observed_h2' or r['plan_sha256']!=expected_hash or r['stage_receipt_sha256']!=stage_hash or r['verifier_sha256']!=sha(Path(__file__)) or r['arithmetic_pass'] is not True:
        raise RuntimeError('EXACT_DIAGNOSTIC_REPORT_IDENTITY_REQUIRED')
    if native['status']!='QUALIFIED_FEASIBILITY_STAGE_COMPLETE_VERIFIED' or native['plan_sha256']!=expected_hash or native['owned_cleanup_verified'] is not True:
        raise RuntimeError('EXACT_COMPLETE_NATIVE_H2_STAGE_REQUIRED')
    native_binding=dict(plan_sha256=expected_hash,admission_sha256=native['admission_sha256'],executor_sha256=sha(executor))
    regular_hashes(native['result_sha256'])
    native_seal=regular(plan['terminal_seal']);native_seal_hash=sha(native_seal)
    require_committed(plan['terminal_pending'],native_seal,native_binding,{str(stage):stage_hash})
    if native_seal_hash!=r['terminal_seal_sha256'] or native['result_sha256'].get(plan['jobs'][0]['result_receipt'])!=r['full_precision_capture_sha256']:
        raise RuntimeError('FROZEN_NATIVE_CAPTURE_OR_TERMINAL_CHANGED')
    spec=importlib.util.spec_from_file_location('_verified_report_prior',executor)
    producer=importlib.util.module_from_spec(spec);spec.loader.exec_module(producer);producer.prior_preprocessing_gate(plan)
    pending,seal=report_terminal_paths(out);regular(seal)
    if r['report_pending_path']!=str(pending) or r['report_terminal_seal_path']!=str(seal):raise RuntimeError('EXACT_REPORT_TERMINAL_ROUTE_REQUIRED')
    pair={str(out):sha(out),str(table):sha(table)}
    binding=dict(plan_sha256=expected_hash,native_stage_receipt_sha256=stage_hash,verifier_sha256=sha(Path(__file__)))
    report_seal_hash=require_committed(pending,seal,binding,pair)
    producer.prior_preprocessing_gate(plan)
    regular_hashes(native['result_sha256']);regular_hashes(pair)
    require_committed(plan['terminal_pending'],native_seal,native_binding,{str(stage):stage_hash})
    if sha(regular(plan_path))!=expected_hash or sha(regular(native_seal))!=native_seal_hash or sha(regular(seal))!=report_seal_hash:
        raise RuntimeError('REPORT_EVIDENCE_CHANGED_DURING_CONSUMPTION')
    return dict(report_pair_sha256=pair,report_terminal_seal_sha256=report_seal_hash,native_terminal_seal_sha256=native_seal_hash)

def verify(plan_path,expected_hash,out):
    plan_path=regular(plan_path)
    if sha(plan_path)!=expected_hash:raise RuntimeError('FROZEN_H2_PLAN_CHANGED')
    plan=json.loads(plan_path.read_text())
    if plan['stage']!='observed_h2' or plan['scope']!='FINNGEN_INSOMNIA_SOURCE_FEASIBILITY_ONLY' or plan['new_rg_commands']!=0 or len(plan['jobs'])!=1:
        raise RuntimeError('SINGLE_QUALIFIED_OBSERVED_H2_SCOPE_REQUIRED')
    job=plan['jobs'][0];prefix=job['out_prefix']
    if job['kind']!='h2' or job['estimates']!=1 or job['inputs']!=[plan['derivative']] or job['output_prefix']!=prefix or job['result_receipt']!=prefix+'.full_precision.json':
        raise RuntimeError('EXACT_SINGLE_CAPTURE_ROUTE_REQUIRED')
    argv=job['ldsc_args']
    if argv.count('--h2')!=1 or argv[argv.index('--h2')+1]!=plan['derivative'] or '--rg' in argv or '--samp-prev' in argv or '--pop-prev' in argv or '--print-delete-vals' not in argv or argv[argv.index('--n-blocks')+1]!='200' or argv[argv.index('--out')+1]!=prefix:
        raise RuntimeError('ONLY_FROZEN_OBSERVED_H2_ARGUMENTS_ADMITTED')
    stage=regular(plan['stage_receipt']);stage_hash=sha(stage);completed=json.loads(stage.read_text())
    executor=regular(plan['executor_path'])
    binding={'plan_sha256':expected_hash,'admission_sha256':completed['admission_sha256'],'executor_sha256':sha(executor)}
    terminal_hash=sha(regular(plan['terminal_seal']))
    require_committed(plan['terminal_pending'],plan['terminal_seal'],binding,{str(stage):stage_hash})
    spec=importlib.util.spec_from_file_location('_qualified_finngen_prior',executor)
    producer=importlib.util.module_from_spec(spec);spec.loader.exec_module(producer)
    producer.prior_preprocessing_gate(plan)
    if completed['status']!='QUALIFIED_FEASIBILITY_STAGE_COMPLETE_VERIFIED' or completed['plan_sha256']!=expected_hash or completed['owned_cleanup_verified'] is not True:
        raise RuntimeError('H2_STAGE_NOT_COMPLETE_AND_REAPED')
    frozen=dict(completed['result_sha256']);regular_hashes(frozen)
    capture=regular(job['result_receipt']);log=prefix+'.log';deletes={prefix+'.delete',prefix+'.part_delete'}
    required={str(capture),log,*deletes,job['worker_receipt']}
    if not required.issubset(frozen):raise RuntimeError('CAPTURE_STOCK_FILES_AND_WORKER_MUST_BELONG_TO_COMPLETED_MAP')
    r=json.loads(capture.read_text())
    if len(r['estimates'])!=1 or len(r['final_intersections'])!=1 or r['plan_sha256']!=expected_hash or r['job_id']!=job['job_id']:
        raise RuntimeError('H2_CAPTURE_IDENTITY_OR_CARDINALITY_CHANGED')
    if r['input_sha256_before']!=plan['input_sha256'] or r['input_sha256_after']!=plan['input_sha256']:
        raise RuntimeError('SOURCE_DERIVATIVE_CHANGED')
    args=r['arguments']
    if args['h2']!=plan['derivative'] or args['rg'] is not None or args['n_blocks']!=200 or args['out']!=prefix or args['samp_prev'] is not None or args['pop_prev'] is not None or args['print_delete_vals'] is not True:
        raise RuntimeError('CAPTURED_OBSERVED_H2_ARGUMENTS_CHANGED')
    if r['stock_log_sha256']!=frozen[log] or set(r['stock_delete_array_sha256'])!=deletes:
        raise RuntimeError('EXACT_STOCK_LOG_AND_TWO_DELETE_HASH_BINDINGS_REQUIRED')
    if any(r['stock_delete_array_sha256'][p]!=frozen[p] for p in deletes):
        raise RuntimeError('CAPTURED_STOCK_DELETE_HASH_DISAGREES_WITH_COMPLETED_MAP')
    worker=json.loads(regular(job['worker_receipt']).read_text())
    if worker['status']!='WORKER_COMPLETE_VERIFIED' or worker['plan_sha256']!=expected_hash or worker['returncode']!=0 or worker['stop_reason'] is not None or worker['owned_cleanup_verified'] is not True or worker['process_group_teardown']['remaining_group_members']:
        raise RuntimeError('EXACT_COMPLETED_REAPED_H2_WORKER_REQUIRED')
    if not {str(capture),log,*deletes}.issubset(worker['output_sha256']) or any(worker['output_sha256'][p]!=frozen.get(p) for p in worker['output_sha256']):
        raise RuntimeError('WORKER_OUTPUT_MAP_DISAGREES_WITH_COMPLETED_STAGE')
    regular_hashes(worker['output_sha256'])
    intersection=r['final_intersections'][0];n=intersection['final_ordered_SNP_count']
    if intersection['input']!=plan['derivative'] or type(n) is not int or n<200 or intersection['n_blocks']!=200:
        raise RuntimeError('ONE_POSITIVE_NATIVE_H2_INTERSECTION_REQUIRED')
    for key in ['final_ordered_SNP_sha256','final_ordered_Z_float64_big_endian_sha256','final_ordered_N_float64_big_endian_sha256']:
        if not re.fullmatch('[0-9a-f]{64}',intersection[key]):raise RuntimeError('NATIVE_INTERSECTION_DIGEST_REQUIRED')
    x=r['estimates'][0]
    if x['input']!=plan['derivative'] or x['n_blocks']!=200 or x['constrain_intercept'] is not False or x['n_annot']!=1:
        raise RuntimeError('SINGLE_INPUT_OR_FREE_INTERCEPT_CHANGED')
    if x.get('n_snp') is not None and x['n_snp']!=n:raise RuntimeError('OBJECT_AND_INTERSECTION_SNP_COUNT_DISAGREE')
    for k in ['tot','tot_se','intercept','intercept_se','mean_chisq','lambda_gc']:
        if type(x[k]) not in (int,float) or not math.isfinite(x[k]):raise RuntimeError('NONFINITE_H2_CAPTURE')
    if x['tot_se']<=0 or x['intercept_se']<=0:raise RuntimeError('NONPOSITIVE_STANDARD_ERROR')
    out=Path(out)
    if out.suffix!='.json' or not out.parent.is_dir() or out.parent.is_symlink() or any(parent.is_symlink() for parent in out.parents):
        raise RuntimeError('ONE_REGULAR_PRIVATE_REPORT_DIRECTORY_REQUIRED')
    pending,seal=report_terminal_paths(out)
    report_binding=dict(plan_sha256=expected_hash,native_stage_receipt_sha256=stage_hash,verifier_sha256=sha(Path(__file__)))
    terminal=TerminalCommit(pending,seal,report_binding)
    checks={};stock_checks={} 
    for suffix,field in [('.delete','tot_delete_values'),('.part_delete','part_delete_values')]:
        native=scalar_deletes(x[field]);stock=stock_scalar_deletes(prefix+suffix)
        passed=all(math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-15) for a,b in zip(native,stock))
        stock_checks[field]=dict(rows=200,fixed_tolerance_pass=passed,stock_file_sha256=frozen[prefix+suffix])
        if not passed:raise RuntimeError('JSON_AND_ACTUAL_STOCK_DELETE_VECTOR_DISAGREE')
    for deletion,field in [('tot_delete_values','tot_se'),('intercept_delete_values','intercept_se')]:
        value=centered_se(scalar_deletes(x[deletion]))
        checks[field]=dict(independent=value,native=x[field],fixed_tolerance_pass=math.isclose(value,x[field],rel_tol=1e-12,abs_tol=1e-15))
    if not all(c['fixed_tolerance_pass'] for c in checks.values()):raise RuntimeError('H2_DELETE_ARITHMETIC_FAILURE')
    z=x['tot']/x['tot_se']
    diagnostic=dict(trait='finngen_R13_F5_INSOMNIA',source='R13_generation_1777989563097164',
        scale='OBSERVED_ONLY_WITH_ASSUMED_CONSTANT_EFFECTIVE_N',h2=x['tot'],h2_SE=x['tot_se'],h2_Z=z,
        h2_CI_lower=x['tot']-1.959963984540054*x['tot_se'],h2_CI_upper=x['tot']+1.959963984540054*x['tot_se'],
        intercept=x['intercept'],intercept_SE=x['intercept_se'],lambda_GC=x['lambda_gc'],mean_chi_square=x['mean_chisq'],
        native_SNP_count=n,assumed_effective_N=plan['assumed_effective_N'],necessary_h2_Z_ge4=z>=4,necessary_intercept_le1_2=x['intercept']<=1.2,
        source_INFO_gate_verified=False,per_variant_N_verified=False,scientific_source_admitted=False,
        fully_independent_two_trait_replication=False,new_rg_estimates=0,
        diagnostic_scope='Source feasibility only; threshold conditions do not resolve INFO/N, ascertainment, overlap, clinical construct or power.')
    result=dict(schema='single_observed_h2_arithmetic_stock_closure_v5',plan_sha256=expected_hash,stage_receipt_sha256=stage_hash,
        full_precision_capture_sha256=frozen[str(capture)],verifier_sha256=sha(Path(__file__)),
        fixed_tolerances=dict(relative=1e-12,absolute=1e-15),arithmetic_checks=checks,stock_delete_checks=stock_checks,
        arithmetic_pass=True,diagnostic=diagnostic,warning_lines=r['warnings_and_errors'],
        intercept_delete_export_scope='Intercept arithmetic uses captured native array; stock LDSC does not export a separate intercept-delete file.',
        all_200_deletions_preserved_in_capture=True,common_cross_trait_boundaries_established=False,original_families_changed=False,
        terminal_seal_sha256=terminal_hash,current_PENDING_or_failure_veto_checked=True,
        report_pending_path=str(pending),report_terminal_seal_path=str(seal),
        report_is_provisional_until_private_PENDING_absent_and_exact_pair_terminal_verified=True)
    def final_identity():
        producer.prior_preprocessing_gate(plan)
        require_committed(plan['terminal_pending'],plan['terminal_seal'],binding,{str(stage):stage_hash})
        if sha(regular(plan['terminal_seal']))!=terminal_hash or sha(regular(plan_path))!=expected_hash:
            raise RuntimeError('TERMINAL_OR_PLAN_CHANGED_AFTER_ARITHMETIC')
        regular_hashes(frozen)
    final_identity()
    table=Path(out).with_suffix('.tsv');buf=io.StringIO(newline='')
    writer=csv.DictWriter(buf,fieldnames=list(diagnostic),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerow(diagnostic)
    payloads={Path(out):json.dumps(result,indent=2,allow_nan=False)+'\n',table:buf.getvalue()}
    expected={str(path):hashlib.sha256(payload.encode()).hexdigest() for path,payload in payloads.items()}
    for path,payload in payloads.items():
        with path.open('x') as f:f.write(payload);f.flush();os.fsync(f.fileno())
    fixed_report_seal_intent=hashlib.sha256((json.dumps(dict(status='REVIEWED_STAGE_TERMINAL_SEAL',
        binding=report_binding,result_receipt_sha256=expected,pending_path=str(pending),pending_sha256=terminal.pending_sha,
        success_requires_absent_PENDING_and_all_failure_addenda=True),indent=2,allow_nan=False)+'\n').encode()).hexdigest()
    def report_identity():
        final_identity();regular_hashes(expected);regular_hashes({str(pending):terminal.pending_sha})
        if seal.exists() or seal.is_symlink():regular_hashes({str(seal):fixed_report_seal_intent})
    if not terminal.commit(expected,identity_gate=report_identity,resource_gate=lambda:None,termination_gate=lambda:None):
        raise RuntimeError('DIAGNOSTIC_REPORT_PAIR_TERMINAL_COMMIT_FAILED_PENDING_VETOES_PROVISIONAL_REPORTS')
    try:print(json.dumps(dict(arithmetic_pass=True,stock_total_and_partition_deletes_agree=True,source_admitted=False,new_rg_estimates=0,report_pair_terminal_complete=True)))
    except BaseException:pass

def main():
    p=argparse.ArgumentParser();p.add_argument('--plan',type=Path,required=True);p.add_argument('--plan-sha256',required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();verify(a.plan,a.plan_sha256,a.out)

if __name__=='__main__':main()

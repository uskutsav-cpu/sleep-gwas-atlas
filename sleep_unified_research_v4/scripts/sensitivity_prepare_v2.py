#!/usr/bin/env python3
"""Additive result-free baseline-gate correction; preserve every sealed v1 file."""
import copy
import datetime
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil

ROOT=Path(__file__).resolve().parents[2]
P=ROOT/'sleep_unified_research_v4'
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/sensitivities')
V1_SHA='cb2f39934851541fe087fab5bac83d332ca32b6b6861d1973da6470b1e06844c'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def main():
    old=SSD/'sensitivity_operational_plan_v1.json'
    if sha(old)!=V1_SHA:raise RuntimeError('SEALED_V1_PLAN_CHANGED')
    plan=copy.deepcopy(json.loads(old.read_text()))
    for path,digest in plan['dependencies_sha256'].items():
        if sha(path)!=digest:raise RuntimeError('SEALED_V1_DEPENDENCY_CHANGED: '+path)
    for path,digest in plan['input_sha256'].items():
        if sha(path)!=digest:raise RuntimeError('SOURCE_CHANGED: '+path)
    evidence={
        'logs/core_native_comparison_receipt_v4.json':'a71b70f18090fd10440e729b1c2df872b4a8e1a2c8208d401db5b9b7b81a81d9',
        'reviews/independent_core_numerical_adjudication_v4.json':'623b6e4e4815286f81a1a52cb98f30ba329bab8438a86d23eed020966302d40d',
        'reviews/independent_core_numerical_adjudication_v4.md':'7fed0ba8c0d0e808a3a02122deedc82247a1ff6ef1020fee5b279ccac83258c7',
        'reviews/independent_whole_core_receipt_v4_2.json':'72732c58a0de3fd87995e748e111d3cf24cc4b157fdde443a59a19b252850fa0',
        'reviews/independent_core_command_binding_receipt_v4.json':'1b1a9b3e93a177dfaeab327d2d68cc0f2bf77a5df4265a12817dee699333ea40',
        'reviews/independent_whole_core_plan_v4_2.json':'4eab44b6c2607fac22daee752d7c81d0990680049c19ec815f15d85bfd008eb0',
        'reviews/independent_whole_core_checker_v4_2.py':'961a5007205e6c2acd30c6875dff10a38d6fe6a0813f7303e99e27a4653eac21',
        'reviews/independent_core_subset_checker_v4_3.py':'bd9ea033b1bc6e4a2347f672e944037541987bc0f747f1eab56e3b983f6f2d70',
        'reviews/independent_core_command_binding_checker_v4.py':'86e8238158a12bb9589e3e984899f089b62d00c1083476c34fd2b05b291e7aae',
        'scripts/33_compare_native_campaign_original_pseudovariance_v1.py':'33ceb0d7e4fc0f3bbb93cb4d04fbf46cf82b24fad942461a2aa8b4e36fa90f44',
    }
    bindings={str(P/path):digest for path,digest in evidence.items()}
    for path,digest in bindings.items():
        if sha(path)!=digest:raise RuntimeError('INDEPENDENT_ADJUDICATION_IDENTITY_DIFFERS: '+path)
    comparison=json.loads((P/'logs/core_native_comparison_receipt_v4.json').read_text())
    if len(comparison['arithmetic_failures'])!=1:raise RuntimeError('CORE_FAILURE_CARDINALITY_NOT_EXACTLY_ONE')
    exact=comparison['arithmetic_failures'][0]
    if exact['identity']!='snoring__bmi' or exact['p'] is not False or any(exact[k] is not True for k in ['hsq1_se','hsq2_se','gencov_se','rg_ratio','rg_jknife','rg_se','z','capture_to_native_summary_pass']):
        raise RuntimeError('CORE_FAILURE_IS_NOT_THE_DECLARED_EXTREME_TAIL_P_RECORD')
    plan['core_precision_adjudication']=dict(
        comparison=str(P/'logs/core_native_comparison_receipt_v4.json'),
        adjudication=str(P/'reviews/independent_core_numerical_adjudication_v4.json'),
        report=str(P/'reviews/independent_core_numerical_adjudication_v4.md'),
        whole_core=str(P/'reviews/independent_whole_core_receipt_v4_2.json'),
        command_binding=str(P/'reviews/independent_core_command_binding_receipt_v4.json'),
        evidence_sha256=bindings,exact_preserved_failure=exact,
        scope='Only preserved snoring__bmi pseudovariance extreme-tail P roundoff; independently centered arithmetic passes unchanged fixed tolerances. No other failure admitted.')
    original=str(ROOT/'sleep_unified_research_v1/scripts/native_ldsc_capture.py')
    relocated=str(Path(plan['baseline_support_package'])/'scripts/native_ldsc_capture.py')
    expected=dict(plan['baseline_dependency_sha256'])
    pinned='bd1b87e403fca2bad4c30b6a662a64f7c3d1a47f6aa9884f5157e370cb95707d'
    if expected.get(original)!=pinned or sha(original)!=pinned or sha(relocated)!=pinned or relocated in expected:
        raise RuntimeError('CAPTURE_RELOCATION_HASH_DIFFERS')
    expected[relocated]=expected.pop(original)
    plan['baseline_relocated_dependency_sha256']=expected
    plan['baseline_dependency_relocation']=dict(original_path=original,relocated_path=relocated,unchanged_sha256=pinned,
        changed_path_keys=1,other_dependency_keys_and_hashes_unchanged=True,source='Original v1 runner execution_dependencies uses relocated PACKAGE only.')
    plan['dependencies_sha256'].update(bindings)
    plan['dependencies_sha256'][relocated]=pinned
    plan['dependencies_sha256'][str(old)]=V1_SHA
    names=['sensitivity_executor_v2.py','sensitivity_prepare_v2.py','verify_sensitivity_corrections_v2.py']
    for name in names:plan['dependencies_sha256'][str(P/'scripts'/name)]=sha(P/'scripts'/name)
    for j in plan['jobs']:
        j['out_prefix']=str(SSD/'fits_v2'/j['job_id']);i=j['ldsc_args'].index('--out');j['ldsc_args'][i+1]=j['out_prefix']
    for j in plan['audit_jobs']:
        j['out_prefix']=str(SSD/'proofs/intersections_v2'/j['audit_id']);i=j['ldsc_args'].index('--out');j['ldsc_args'][i+1]=j['out_prefix']
    for directory in ['fits_v2','logs_v2','receipts_v2','proofs/intersections_v2']:(SSD/directory).mkdir(parents=True,exist_ok=True)
    spec=importlib.util.spec_from_file_location('_result_free_v2_gate_check',P/'scripts/sensitivity_executor_v2.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.relocated_baseline_dependencies(plan);mod.admit_exact_core_adjudication(plan,comparison)
    plan['schema']='frozen_estimator_sensitivity_operational_plan_v2'
    plan['prepared_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    plan['superseded_operational_plan']=str(old);plan['superseded_operational_plan_sha256']=V1_SHA
    plan['executor']=str(P/'scripts/sensitivity_executor_v2.py')
    plan['operational_corrections']=['Exact one-key hash-bound baseline capture relocation','Exactly one independently adjudicated preserved core P-arithmetic failure','Unconditional owned teardown before guarded final queries and primary/fallback failure receipts']
    plan['status']='PREPARED_V2_ONLY_NO_AUDITS_OR_FITS_LAUNCHED'
    plan['resource_preflight']=dict(internal_free_bytes=shutil.disk_usage('/System/Volumes/Data').free,SSD_free_bytes=shutil.disk_usage(SSD).free)
    if plan['resource_preflight']['internal_free_bytes']<3*(1<<30) or plan['resource_preflight']['SSD_free_bytes']<5*(1<<30):raise RuntimeError('ORIGINAL_RESOURCE_FLOORS_FAIL')
    out=SSD/'sensitivity_operational_plan_v2.json'
    with out.open('x') as f:json.dump(plan,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(plan=str(out),sha256=sha(out),dependencies=len(plan['dependencies_sha256']),native_commands_prepared=26,new_estimates_prepared=62,launched_workers=0),indent=2))


if __name__=='__main__':main()

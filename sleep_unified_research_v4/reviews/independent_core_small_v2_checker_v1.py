#!/usr/bin/env python3
"""Narrow, metadata-only review of the frozen core-small v2 command amendment.
No production GWAS/reference bodies, complete runtime trees, workers or fits.
"""
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

sys.dont_write_bytecode = True
P = Path(__file__).resolve().parents[1]
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OLD = SSD / 'core_pipeline/original_small_input_replay_v1'
NEW = SSD / 'core_pipeline/original_small_input_replay_v2'
PLAN = NEW / 'core_original_small_replay_plan_v2.json'
PROOF = P / 'reviews/independent_core_small_v2_receipt_v1.json'
OWN = SSD / 'tmp/independent_core_small_v2_metadata_controls_v1'
EXPECTED_PLAN = '9bc40795f552ea203a1a307cdef1e6c5c5f43442407be2e1a559a22325389463'
CHECKS = []
BINDINGS = {}
T0 = time.monotonic()

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()

def bind(p, expected=None):
    got=sha(p)
    assert expected is None or got==expected,(str(p),got,expected)
    BINDINGS[str(Path(p))]=got
    return got

def check(name, good, **detail):
    CHECKS.append(dict(name=name,pass_control=bool(good),**detail))
    assert good,CHECKS[-1]

def load(p,expected=None):
    bind(p,expected)
    return json.loads(Path(p).read_text())

def main():
    assert not PROOF.exists() and not OWN.exists()
    OWN.mkdir(parents=True)
    plan=load(PLAN,EXPECTED_PLAN)
    old=load(OLD/'core_original_small_replay_plan_v1.json','4b63e99da34f7ddcb05ab7109882c5caefea86d913c52ba417cf2e414e9d8b76')
    previous=load(P/'reviews/core_original_small_replay_independent_review_v1.json')
    seal=load(P/'reviews/core_original_small_replay_independent_review_seal_v1.json')
    for path,digest in seal['file_sha256'].items():bind(path,digest)
    controls=load(previous['controls_receipt_path'],previous['controls_receipt_sha256'])
    check('inherited_41_controls_are_passed_metadata_stubs', controls['all_control_expectations_pass']
          and controls['control_count']==41 and controls['actual_workers']==0 and controls['fits']==0)
    terminalseal=load(P/'reviews/terminal_commit_common_review_seal_v2.json')
    for path,digest in terminalseal['file_sha256'].items():bind(path,digest)
    compat=load(P/'reviews/independent_historical_munge_text_receipt_v1.json',
                'a8e23d49b8bde759598f761757d8c5a242ecc5d25cb5b3a2b2de16e439de5de6')
    bind(P/'reviews/independent_historical_munge_text_review_v1.md',
         '1168ebb394e1579082d9283e4df2a3e026c4ecea5bc7ff59ceaa4051d0e354b1')
    check('compatibility_qualified_not_full_replay_or_binary_attestation',
          compat['status']=='QUALIFIED_PRELAUNCH_METHODS_PASS'
          and compat['passed_control_count']==43 and not compat['scientific_full_replay_certified']
          and not compat['historical_per_trait_binary_attestation'] and not compat['execution_admission'])
    provenance=load(plan['historical_munge_compatibility_provenance'],
                    plan['historical_munge_compatibility_provenance_sha256'])
    candidate=provenance['candidate_executable']
    oldexe=provenance['pinned_unmodified_munger']
    check('exact_candidate_hash_and_qualifier_bound',
          plan['dependencies_sha256'][candidate]==provenance['candidate_executable_sha256']
          and plan['historical_munge_qualification']==provenance['qualification']
          and 'no historical binary identity claim' in plan['execution_preconditions'])
    # References are named explicitly and are not opened, hashed or decompressed.
    excluded=[]
    for path,digest in plan['dependencies_sha256'].items():
        f=Path(path)
        if f.name in ['hm3_grch37_variant_map.tsv.gz','hg38ToHg19.over.chain.gz','w_hm3.snplist']:
            excluded.append(dict(path=path,expected_sha256=digest,bytes=f.stat().st_size));continue
        assert f.suffix not in ['.gz','.bz2','.bgz'],('unexpected body dependency',path)
        bind(f,digest)
    check('exact_dependency_metadata_set',len(plan['dependencies_sha256'])==167 and len(excluded)==3,
          verified=164,unread_reference_bodies=excluded)
    names=['ms','asthma','bmi','t2d','ldl','hdl','triglycerides','cad','telomere_length','melanoma']
    check('original_ten_members_and_35_remaining', [m['trait_id'] for m in plan['members']]==names
          and plan['member_count']==10 and plan['remaining_large_input_traits']==35
          and plan['original_core_family_size']==45)
    check('original_design_rows_unchanged', [m['original_design'] for m in plan['members']]==
          [m['original_design'] for m in old['members']])
    def relocate(value):
        if value==oldexe:return candidate
        return str(NEW)+value[len(str(OLD)):] if value.startswith(str(OLD)+'/') else value
    for original,current in zip(old['members'],plan['members']):
        for key in ['prefilter_command','harmonize_command','munge_command']:
            expected=None if original[key] is None else [relocate(x) for x in original[key]]
            check('exact_command_amendment_'+current['trait_id']+'_'+key,current[key]==expected)
        mu=current['munge_command']
        assert mu[2]==candidate and mu[mu.index('--chunksize')+1]=='500000'
        for key in ['harmonized','harmonization_qc','munged','munged_prefix','source_gate_receipt','comparison_receipt']:
            assert current[key]==relocate(original[key]),key
    stable=['baseline_support_package','baseline_execution_plan_sha256','baseline_dependency_sha256',
            'baseline_relocated_dependency_sha256','baseline_input_sha256','original_190_jobs',
            'core_precision_adjudication','python','ldsc_dir','reference_prefix','environment','guard',
            'runtime_qualification','qualified_runtime_receipt','qualified_runtime_receipt_sha256',
            'harmonization_python','scientific_membership_or_threshold_changes','automatic_retry',
            'estimator_calls','archived_input_sha256']
    check('baseline190_options_inputs_estimators_runtime_and_guards_unchanged',
          all(plan[k]==old[k] for k in stable),unchanged_keys=stable)
    check('exact_46_command_scope',sum(m['prefilter_command'] is not None for m in plan['members'])==6
          and len(plan['original_190_jobs'])==190 and plan['estimator_calls']==0,
          command_counts={'source_validation':10,'prefilter':6,'harmonize':10,'munge':10,'comparison':10,'total':46})
    executor=P/'scripts/71_run_core_original_small_replay_v2.py'
    prior=P/'scripts/71_run_core_original_small_replay_v1.py'
    txt=executor.read_text()
    normalized=txt
    for before,after in [('original_small_input_replay_v2','original_small_input_replay_v1'),
                         ('core_pending_v2','core_pending_v1'),('core_terminal_seal_v2','core_terminal_seal_v1'),
                         ('72_validate_core_original_small_replay_v2','72_validate_core_original_small_replay_v1'),
                         ('core_original_small_execution_receipt_v2','core_original_small_execution_receipt_v1')]:
        normalized=normalized.replace(before,after)
    check('controller_whole_text_identity_except_five_private_path_literals',normalized==prior.read_text()
          and ast.dump(ast.parse(normalized),include_attributes=False)==ast.dump(ast.parse(prior.read_text()),include_attributes=False))
    check('validator_whole_byte_identity',sha(P/'scripts/72_validate_core_original_small_replay_v2.py')==
          sha(P/'scripts/72_validate_core_original_small_replay_v1.py'))
    failure=load(plan['preserved_failed_v1_receipt'],plan['preserved_failed_v1_receipt_sha256'])
    check('failed_v1_immutable_pending_and_worker_preserved',failure['status']=='FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
          and failure['owned_cleanup_verified'] and (OLD/'core_pending_v1.json').is_file()
          and str(OLD/'core_pending_v1.json') in plan['dependencies_sha256']
          and str(OLD/'receipts_v4/ms__munge.worker.json') in plan['dependencies_sha256']
          and str(OLD/'logs_v4/ms__munge.worker.stdout.log') in plan['dependencies_sha256'])
    ledger=load(plan['global_resource_ledger_path'],plan['global_resource_ledger_sha256'])
    check('ledger3_counts_old_failures_and_16GiB_core',sum(ledger['component_bytes'].values())==ledger['reserved_total_bytes']
          and ledger['reserved_total_bytes']+ledger['unallocated_margin_bytes']==ledger['ceiling_bytes']==300*2**30
          and ledger['component_bytes']['proposed_core_derivative_cap_bytes']==16*2**30
          and ledger['core16GiB_cap_includes_final_outputs_spools_indexes_prefilters_transient_and_failed_attempts'])
    check('post_persistence_final_identity_resource_termination_gate_inherited',
          'identity_gate=final_identity' in txt and "POST_PERSISTENCE_CORE_TERMINAL" in txt
          and "TerminalCommit(OUT/'core_pending_v2.json',OUT/'core_terminal_seal_v2.json'" in txt
          and plan['terminal_contract']==old['terminal_contract'])
    # Invoke only the real admission gate with invented metadata; rejection is
    # before check_bindings/runtime/body reads or any worker construction.
    sys.path.insert(0,str(P/'scripts'))
    spec=importlib.util.spec_from_file_location('_independent_core_small_v2_admission',executor)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    for name,admission in [
      ('not_admitted',dict(execution_admitted=False)),
      ('missing_review',dict(execution_admitted=True,plan_sha256=EXPECTED_PLAN,executor_sha256=sha(executor),independent_review_sha256={})),
      ('review_hash_changed',dict(execution_admitted=True,plan_sha256=EXPECTED_PLAN,executor_sha256=sha(executor),
                                independent_review_sha256={str(P/'reviews/independent_historical_munge_text_receipt_v1.json'):'0'*64}))]:
        path=OWN/(name+'.json')
        with path.open('x') as f:json.dump(admission,f)
        args=SimpleNamespace(plan=PLAN,plan_sha=EXPECTED_PLAN,admission=path,admission_sha=sha(path))
        error=None
        try:mod.admission_gate(plan,args)
        except RuntimeError as e:error=str(e)
        check('real_admission_rejection_'+name,error is not None,reason=error)
        bind(path)
    for path,digest in BINDINGS.items():assert sha(path)==digest,('after_review_mutation',path)
    result={'schema':'independent_narrow_core_small_v2_prelaunch_review_v1',
            'verdict':'QUALIFIED_NARROW_PRELAUNCH_PASS', 'plan_path':str(PLAN),'plan_sha256':EXPECTED_PLAN,
            'checks':CHECKS,'check_count':len(CHECKS),'consumed_metadata_code_before_after_sha256':BINDINGS,
            'metadata_dependencies_verified':164,'reference_body_dependencies_unread':excluded,
            'inherited_controller_control_count':41,'inherited_controls_rerun':False,
            'scientific_arguments_order_and_values_preserved':True,'historical_munge_header_text_qualification_required':True,
            'execution_admission_granted':False,'source_or_reference_body_reads':0,'real_data_decompressions':0,
            'actual_workers':0,'fits':0,'material_new_blockers':[],
            'not_certified':['future actual full ordered content and QC replay','fresh runtime-tree/source/reference identities',
                             'historical per-trait binaries','actual production peak RSS','future downstream terminal consumer compliance'],
            'required_root_admission':['bind this exact plan and executor and this narrow MD/JSON review',
                                       'bind exact historical-munger compatibility MD/JSON review and retain its qualifier',
                                       'require complete original190 numerical gate and fresh exact runtime/source checks',
                                       'success requires exact terminal seal/current receipt hashes/absent PENDING and failure addenda'],
            'elapsed_seconds':time.monotonic()-T0}
    with PROOF.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps({'verdict':result['verdict'],'checks':len(CHECKS),'receipt_sha256':sha(PROOF),'elapsed':result['elapsed_seconds']}))

if __name__=='__main__':main()

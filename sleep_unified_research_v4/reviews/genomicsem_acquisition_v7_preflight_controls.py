#!/usr/bin/env python3
"""Narrow v7 metadata/origin/cardinality controls; never hashes a raw body."""
import ast
import copy
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import sys
import tempfile
import time
from types import SimpleNamespace
from unittest.mock import patch

P=Path(__file__).resolve().parents[1]
S=P/'scripts'
PLAN=P/'manifests/extension_raw_acquisition_plan_v4_7.json'
OLD_PLAN=P/'manifests/extension_raw_acquisition_plan_v4_6.json'
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
EXECUTOR=S/'52_acquire_extension_raw_sources_v7.py'
PINS={str(PLAN):'0958a8383ada38c25eceb6451a7d074527753fae5dffce278d02e8bb2ba1345c',
      str(EXECUTOR):'f4f5e7b8c5c1f4942ff65a95737b525c1df9f9f53633c1fa269a1c774976de12',
      str(OLD_PLAN):'393af7ffc3231ceddf829ae669558001a7411bb4a7917b950a713c41546237fe',
      str(P/'manifests/global_SSD_resource_reservation_v4_3.json'):'6bd5c7cec9d9880ae630e5bfd431b548bd4cc2f3837d3ce12935fda326a9c590'}


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()


def module(path,name):
    sp=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m


def funcs(path):
    return {n.name:n for n in ast.parse(path.read_text()).body if isinstance(n,ast.FunctionDef)}


def run():
    checks=[]
    def check(value,label):
        if not value:raise AssertionError(label)
        checks.append({'case':label,'pass':True})
    for p,h in PINS.items():check(sha(p)==h,'pin '+p)
    plan=json.loads(PLAN.read_text());old=json.loads(OLD_PLAN.read_text())
    forbidden={m['body_path'] for x in [plan,old] for m in x['members']}
    forbidden.update(p+'.partial' for p in list(forbidden));forbidden.add(old['explicit_resume']['original_partial_path']);forbidden.add(plan['explicit_resume']['original_partial_path'])
    before={}
    for p,h in plan['bound_sources'].items():
        check(p not in forbidden and not Path(p).name.endswith(('.gz','.bgz','.bed','.bim','.fam')),'non_body_binding '+p)
        check(not Path(p).is_symlink() and Path(p).is_file() and sha(p)==h,'bound_metadata '+p)
        before[p]={'sha256':h,'bytes':Path(p).stat().st_size}
    copy_plan=SSD/'manifests'/PLAN.name
    check(sha(copy_plan)==PINS[str(PLAN)],'exact_SSD_plan_copy')
    check(len(plan['members'])==len(old['members'])==100,'exact100_members')
    for a,b in zip(old['members'],plan['members']):
        expected=dict(a)
        if a['index']>3:expected['body_path']=str(SSD/'extension_raw_replay_v7/raw'/a['filename'])
        check(expected==b,'ordered_original_source_identity '+b['extension_trait_id'])
    check(set(plan['reused_checkpoints'])=={'1','2','3'} and plan['network_sources_to_transfer']==97,'exact3_reused97_new')
    check(all(plan['bound_sources'].get(p)==h for p,h in old['bound_sources'].items()),'all137_v6_bindings_retained')
    for key in ['compressed_network_bytes','SSD_reservation_bytes','largest_body_bytes','internal_floor_bytes','ssd_floor_bytes','maximum_owned_transfer_rss_bytes','per_body_seconds_limit','family_seconds_limit','transfer_worker_count','automatic_retry','retained_raw_bodies','scientific_membership_or_threshold_changes','curl_path','curl_version','monitor_path','monitor_output_limit_bytes','runtime_poll_seconds','exclusive_family_lock_path','terminal_protocol','all_failed_provisional_primary_copies_must_be_rejected']:
        check(plan[key]==old[key],'unchanged_operational_scientific_policy '+key)
    check(plan['automatic_retry'] is False and plan['transfer_worker_count']==1,'one_transfer_no_automatic_retry')
    r=plan['explicit_resume'];check(r['source_index']==4 and r['source_trait']==plan['members'][3]['extension_trait_id'],'source4_exact_trait')
    check(r['prefix_bytes']==1126899820 and r['prefix_sha256']=='77a3c9be4197fe2f45b529d271effd40f140e117a652079a66870b4fc4743e1d','source4_exact_preserved_prefix')
    check(r['original_partial_path']==old['members'][3]['body_path']+'.partial','original_prefix_in_failed_v6_namespace')
    check(plan['new_network_bytes']==sum(m['expected_bytes'] for m in plan['members'][3:])-r['prefix_bytes']==219656952268,'exact_remaining_network_bytes')
    check(plan['additional_preserved_prefix_bytes']==old['additional_preserved_prefix_bytes']+r['prefix_bytes']==3532333164,'all_three_extra_prefix_reservations')
    failedpath=P/'logs/extension_raw_acquisition_family_receipt_v4_6.json';failed=json.loads(failedpath.read_text());sourcefail=json.loads(Path(r['prior_failed_receipt']).read_text())
    check(sha(failedpath)==plan['prior_v6_failed_family_sha256'] and failed['status']=='FAMILY_STOP_PRESERVED_REQUIRES_EXPLICIT_AUDIT' and failed['completed_source_count']==3,'v6_failed_family_remains_failed')
    check(sourcefail['returncode']==56 and sourcefail['status']=='FAILED_PRESERVED_NO_AUTOMATIC_RETRY' and sourcefail['stop_reason']=='AssertionError: CURL_NONZERO_EXIT: 56','actual_failed_source4_curl56_metadata')
    check(sourcefail['teardown']['teardown_verified'] and sourcefail['teardown']['remaining_group_members']==[],'failed_source4_recorded_cleanup_empty')
    check(sourcefail['retained_partial_sha256']==r['prefix_sha256'] and sourcefail['retained_partial_bytes']==r['prefix_bytes'],'failed_source4_prefix_receipt_matches')
    check(Path(old['pending_path']).is_file() and not Path(old['terminal_seal_path']).exists(),'v6_pending_retained_no_terminal_credit')
    successful={x['path']:x['sha256'] for x in failed['source_receipts']}
    gatepath=Path(plan['reused_checkpoint_gate_plan']);gateplan=json.loads(gatepath.read_text())
    common=module(S/'extension_replay_common_v3.py','_v7_origin_blueprint')
    accepted=[]
    for member in plan['members'][:3]:
        fixed=plan['reused_checkpoints'][str(member['index'])];receipt=json.loads(Path(fixed['receipt_path']).read_text())
        check(sha(fixed['receipt_path'])==fixed['receipt_sha256']==successful[fixed['receipt_path']],'reused_receipt_exact_hash_and_failed_family_map '+str(member['index']))
        check(receipt['member']==member and receipt['plan_sha256']==fixed['original_acquisition_plan_sha256'],'reused_exact_member_and_origin '+str(member['index']))
        expected=old['reused_checkpoints'].get(str(member['index']))
        check(fixed==expected if expected else fixed['original_acquisition_plan_sha256']==PINS[str(OLD_PLAN)],'preserved_old_origin_or_exact_v6_origin '+str(member['index']))
        oldmember=next(m for m in gateplan['members'] if m['index']==member['index'])
        # Exact gate reads only tiny receipt/header bytes and stats body size;
        # unlike reused_checkpoint_gate, it never hashes the GWAS body.
        check(common.acquisition_receipt_gate(gateplan,oldmember,fixed['receipt_sha256'])==receipt,'actual_three_origin_receipt_header_stat_gate '+str(member['index']))
        accepted.append({'index':member['index'],'receipt_path':fixed['receipt_path'],'sha256':fixed['receipt_sha256'],'origin_plan_sha256':fixed['original_acquisition_plan_sha256'],'body_hashes_independently_reread':False})
    ledgerpath=Path(plan['global_resource_ledger_path']);ledger=json.loads(ledgerpath.read_text());priorledger=json.loads(Path(old['global_resource_ledger_path']).read_text())
    check(sha(ledgerpath)==plan['global_resource_ledger_sha256']==PINS[str(ledgerpath)],'new_global_ledger_identity')
    check(sum(ledger['component_bytes'].values())==ledger['reserved_total_bytes']==319645174442,'new_ledger_components_sum')
    check(ledger['reserved_total_bytes']<ledger['ceiling_bytes']==322122547200 and ledger['unallocated_margin_bytes']==2477372758,'new_ledger_under300GiB')
    check(ledger['reserved_total_bytes']-priorledger['reserved_total_bytes']==r['prefix_bytes'],'ledger_increase_only_preserved_source4_prefix')
    a,b=funcs(S/'52_acquire_extension_raw_sources_v6.py'),funcs(EXECUTOR)
    changed={'prepare','reused_checkpoint_gate','acquire','finalize_family','execute'}
    for name in a:
        if name not in changed:check(ast.dump(a[name],include_attributes=False)==ast.dump(b[name],include_attributes=False),'unchanged_v6_function '+name)
    class ProvenanceOnly(ast.NodeTransformer):
        def visit_Constant(self,node):
            if isinstance(node.value,str):node.value=node.value.replace('extension_raw_acquisition_v4_7','extension_raw_acquisition_v4_6')
            return node
    check(ast.dump(a['acquire'],include_attributes=False)==ast.dump(ProvenanceOnly().visit(copy.deepcopy(b['acquire'])),include_attributes=False),'acquire_AST_identical_except_provenance_version_path')
    # Only the changed3/97 cardinality logic is re-exercised. All persistence,
    # quarantine and cleanup faults remain covered by the exact prior12+8 tests.
    mod=module(EXECUTOR,'_v7_cardinality_control');observations=[]
    with tempfile.TemporaryDirectory(prefix='raw-v7-count-review-') as td:
        for good in [True,False]:
            folder=Path(td)/('good' if good else 'wrong_counts');folder.mkdir();package=folder/'package';(package/'logs').mkdir(parents=True)
            pending=folder/'pending.json';terminal=folder/'terminal.json';mod.write_new(pending,{'synthetic_pending':True})
            source_receipts=[]
            for j in range(100):
                path=folder/('receipt_%03d.json'%j);mod.write_new(path,{'synthetic_source':j});source_receipts.append({'path':str(path),'sha256':sha(path)})
            fixture={'reused_checkpoints':{'1':{},'2':{},'3':{}},'terminal_protocol':'SYNTHETIC_EXACT_PROTOCOL','terminal_seal_path':str(terminal),'pending_path':str(pending),'internal_floor_bytes':3*1024**3,'ssd_floor_bytes':5*1024**3,'family_seconds_limit':1000}
            master={'status':'ALL100_EXACT_SOURCE_BODIES_ACQUIRED','source_receipts':source_receipts,'reused_exact_source_count':3 if good else 2,'new_exact_source_count':97 if good else 98}
            monitor=SimpleNamespace(snapshot=lambda:{'internal_free_bytes':4*1024**3,'ssd_free_bytes':6*1024**3})
            with patch.object(mod,'PACKAGE',package),patch.object(mod,'FOLDER',folder),patch.object(mod,'assert_bindings',lambda *_:None),patch.object(mod,'physical_mount',lambda:None),patch.object(mod,'global_namespace_gate',lambda *_:0),patch('subprocess.Popen',side_effect=RuntimeError('FORBIDDEN_WORKER_LAUNCH')):
                committed=mod.finalize_family(fixture,'SYNTHETIC_EXACT_PLAN',master,monitor,time.monotonic(),pending,sha(pending))
            copies=[package/'logs/extension_raw_acquisition_family_receipt_v4_7.json',folder/'extension_raw_acquisition_family_receipt_v4_7.json']
            if good:
                s=json.loads(terminal.read_text());check(committed and not pending.exists() and s['primary_receipt_sha256']=={str(p):sha(p) for p in copies} and len({sha(p) for p in copies})==1,'healthy3_reused97_new_finalizer_two_exact_copies')
            else:check(not committed and pending.exists() and all(Path(str(p)+'.failure.json').exists() for p in [*copies,terminal]),'wrong2_reused98_new_counts_reject_pending_retained')
            observations.append({'good3_97':good,'committed':committed,'pending_retained':pending.exists(),'real_workers':0,'scope':'Changed cardinality only; identity/mount/namespace callbacks labelled fixture mocks.'})
    for p,meta in before.items():check(sha(p)==meta['sha256'] and Path(p).stat().st_size==meta['bytes'],'bound_metadata_after '+p)
    for p,h in PINS.items():check(sha(p)==h,'pin_after '+p)
    check(sha(copy_plan)==PINS[str(PLAN)],'SSD_plan_copy_unchanged_after')
    return {'recorded_utc':datetime.now(timezone.utc).isoformat(),'status':'NARROW_V7_PREPARED_METADATA_ORIGIN_AND_CARDINALITY_PASS','prepared_plan_sha256':PINS[str(PLAN)],'executor_sha256':PINS[str(EXECUTOR)],'helper_sha256':sha(__file__),'pinned_sha256':PINS,'metadata_bindings_verified_before_and_after':before,'metadata_binding_count':len(before),'metadata_bytes_per_pass':sum(x['bytes'] for x in before.values()),'check_count':len(checks),'checks':checks,'accepted_three_origin_metadata_proofs':accepted,'new_cardinality_controls':observations,'inherited_v6_12_terminal_and8_cleanup_controls_not_repeated':True,'no_raw_reference_body_reads':True,'no_gzip_EOF_or_pipeline_claim':True,'raw_body_whole_hashes_are_root_preparation_claim_not_this_review':True,'real_subprocesses_launched':0,'production_mutex_acquired':False,'estimator_calls':0,'execution_admitted':False,'prior_v6_failure_not_reclassified':True,'self_peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}


if __name__=='__main__':
    out=P/'reviews/genomicsem_acquisition_v7_controls_receipt.json'
    if out.exists():raise RuntimeError('DISTINCT_REVIEW_OUTPUT_REQUIRED')
    result=run()
    with out.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','metadata_binding_count','check_count','self_peak_RSS_bytes']}))

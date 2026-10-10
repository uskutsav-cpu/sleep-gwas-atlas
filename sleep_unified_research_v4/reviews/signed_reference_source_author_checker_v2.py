"""Metadata-only immutable prospective v2 closure; no actual bodies/locks/jobs."""
import ast
import datetime
import hashlib
import json
import os
from pathlib import Path

P=Path(__file__).resolve().parents[1];R=P/'reviews';S=P/'scripts'
PLAN=P/'manifests/signed_reference_source_execution_plan_v2.json'
PLAN_SHA='49b8fed7a8af638362dbbededbb9a01e6877f9b727655a3b7327fdfd07c0ccee'


def sha(p):
    p=Path(p);assert p.is_file() and not p.is_symlink(),str(p)
    return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p,h):
    assert sha(p)==h;value=json.loads(Path(p).read_text());assert sha(p)==h;return value
def write(p,value):
    b=(json.dumps(value,sort_keys=True,indent=2)+'\n').encode();h=hashlib.sha256(b).hexdigest()
    with Path(p).open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
    assert sha(p)==h;return h
def funcs(p):
    return {n.name:ast.dump(n,include_attributes=False) for n in ast.parse(Path(p).read_text()).body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}


def main():
    plan=read(PLAN,PLAN_SHA);deferred=plan['preparation_deferred_dependency_sha256']
    assert len(deferred)==26
    fresh={p:h for p,h in plan['dependency_sha256'].items() if p not in deferred};assert len(fresh)==75
    before={p:sha(p) for p in fresh};assert before==fresh
    old_path=P/'manifests/signed_reference_source_execution_plan_v1.json'
    old=read(old_path,'6d73f45ee3f73f7db7c0cc727d9b2d3e1af30d3ee27391cf1154f1f8ab8482b4')
    same=['candidate','archive_members','reference_only_QC','independent_controls','finite_reference_truth_contract',
          'preserved_full_protocol','expected_master_reference_rows','reference_files_numeric_order','coordinate_map_path',
          'authoritative_panel_path','reference_body_sha256','runtime_A_environment_root','runtime_A_profile_sha256',
          'runtime_B_environment_root','runtime_B_profile_sha256','runtime_B_required_versions','maintained_decoderB_path',
          'maintained_decoderB_sha256','decoderB_contract','source_build_qualification','private_data_policy','terminal_contract']
    for k in same:assert plan[k]==old[k],k
    for k,v in old['resource_contract'].items():
        if k!='source_namespace':assert plan['resource_contract'][k]==v,k
    assert deferred==old['preparation_deferred_dependency_sha256']
    code={str(S/(n+'_v2.py')):fresh[str(S/(n+'_v2.py'))] for n in ['signed_reference_source_prepare',
      'signed_reference_source_common','signed_reference_source_controller','signed_reference_source_worker',
      'signed_reference_decode','signed_reference_ldsc_decoderB','signed_reference_compact']}
    for p in code:ast.parse(Path(p).read_text(),filename=p)
    unchanged={}
    for stem,names in [('signed_reference_decode',list(funcs(S/'signed_reference_decode_v1.py'))),
      ('signed_reference_ldsc_decoderB',list(funcs(S/'signed_reference_ldsc_decoderB_v1.py'))),
      ('signed_reference_source_common',['record_output','load_frozen','save','binary_new','streamed_new','full_tree_bytes','Guard','check_dependencies','identity_gate','shared_raw_family_lock','verify_outputs']),
      ('signed_reference_source_controller',['monitor_helper','group_exists','cleanup','runtime_profile_gate','monitor_worker']),
      ('signed_reference_source_worker',['line','consume_jsonl','acquire','configure_db','source_format','read_dose'])]:
        a=funcs(S/(stem+'_v1.py'));b=funcs(S/(stem+'_v2.py'))
        for n in names:assert a[n]==b[n],(stem,n)
        unchanged[stem]=names
    assert plan['preparation_only'] is True and plan['execution_admitted'] is False
    assert plan['prior_v1_rejection_seal_sha256']=='4b024a5f1583b5d7db08f1b20b0cad6f211946a83dde9788d93af6f501b78e8b'
    cap=plan['compact_representation']['capacity']
    assert cap['maximum_planned_QC_bytes']==416822968 and cap['unchanged_QC_cap_bytes']==536870912 and cap['structural_margin_bytes']==120047944
    latest=R/'signed_reference_source_delta_controls_v2_3/delta_control_receipt.json'
    ctl=read(latest,'cc0866656a1a5e4c0ce08a7a0e5976a9d8839ab592396c6396e61b87b90bcc22')
    assert ctl['case_count']==23 and all(c['status']=='PASS' for c in ctl['cases'])
    inventory={'regular_file_sha256':{},'symlink_literal':{},'physical_directories':[],
      'scope':'ALL INVENTED FIXTURES; NO REAL SOURCE/FAM/SNP/INDIVIDUAL EXPORTS'}
    for name in ['signed_reference_source_delta_controls_v2','signed_reference_source_delta_controls_v2_2','signed_reference_source_delta_controls_v2_3']:
        for base,dirs,files in os.walk(R/name,followlinks=False):
            inventory['physical_directories'].append(base)
            for n in dirs+files:
                p=Path(base)/n
                if p.is_symlink():inventory['symlink_literal'][str(p)]=os.readlink(p)
                elif p.is_file():inventory['regular_file_sha256'][str(p)]=sha(p)
    inventory['physical_directories'].sort()
    inv_sha=write(R/'signed_reference_source_delta_inventory_v2.json',inventory)
    after={p:sha(p) for p in fresh};assert before==after==fresh and sha(PLAN)==PLAN_SHA
    receipt=dict(schema='PROSPECTIVE_SOURCE_V2_AUTHOR_PREPARATION_CLOSURE',prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
      verdict='AUTHOR_FOCUSED_DELTA_PASS_CONDITIONAL_CAPACITY_NO_SOURCE_EXECUTION_ADMISSION',plan_sha256=PLAN_SHA,
      new_seven_code_sha256=code,fresh_nonbody_dependency_sha256_before=before,fresh_nonbody_dependency_sha256_after=after,
      inherited26_deferred_identities_unread_and_not_currently_hashed=deferred,unchanged_scientific_plan_fields=same,
      unchanged_AST_scientific_runtime_monitor_cleanup_identity_functions=unchanged,
      current23_delta_control_receipt_sha256=sha(latest),completed42_v1_groups_inherited_not_repeated=True,
      v1_formal_rejection11_and_author53_maps_and_root_consumption_bound_in_actual_plan=True,
      all_four_v1_logical_streams_exact_and_frequency_MAF_bits_checked=True,
      invented_fixture_inventory_sha256=inv_sha,invented_regular_file_count=len(inventory['regular_file_sha256']),
      conservative_completed_retention_capacity=cap,
      capacity_scope='CONDITIONAL_COMPLETE_RETENTION_ALLOCATION_AND_FAIL_CLOSED_GATES; NOT_ACTUAL_SOURCE_FIT_OR_INSTANTANEOUS_FAILED_PARTIAL_PEAK_CERTIFICATION',
      ancillary_controls_metadata_transport_failures_observed_caps_guarded_not_claimed_unconditional_peak_bounds=True,
      IDpool_SQLite_limits_fail_whole_stage_before_genotype_QC_no_dropping_truncation_or_eligibility_change=True,
      actual_eligible_J_count_observed=False,actual_source_archive_genotype_reference_body_operations=0,
      actual_transfer_production_workers_or_mutex_operations=0,S_simulations_run=0,C_calibration_credit=0,estimator_calls=0,
      unchanged_2GiB_RSS_2hour_3GiB_source_512MiB_QC_300GiB_global_guards=True,
      separate_independent_whole_actual_v2_review_and_root_admission_required=True,root_source_execution_admission=False,
      preserved_v1_hypothetical_presentation_rounding='V1_JSON_CEILING832143944;V1_MD_TEXT_TRUNCATION832143943;BOTH_IMMUTABLE_NO_SCIENTIFIC_EFFECT')
    h=write(R/'signed_reference_source_author_metadata_receipt_v2.json',receipt)
    print(json.dumps(dict(status=receipt['verdict'],receipt_sha256=h,new_seven_codes=7,fresh_nonbody=75,deferred=26,delta_groups=23)))


if __name__=='__main__':main()

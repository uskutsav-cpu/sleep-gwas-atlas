"""Metadata/code-only fixed preparation. Never launches/downloads/hashes bodies."""
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path

from signed_reference_source_common_v2 import (P,PLAN,ADMISSION,DESIGN,DESIGN_SHA,
    DESIGN_SEAL,DESIGN_SEAL_SHA,LEDGER,LEDGER_SHA,load_frozen,save,regular_sha)
from signed_reference_compact_v2 import CAPACITY, SCHEMA, SCHEMA_SHA256


def prepare():
    if PLAN.exists() or ADMISSION.exists():raise RuntimeError('FROZEN_SOURCE_PREPARATION_NO_OVERWRITE')
    design=load_frozen(DESIGN,DESIGN_SHA);sealed=load_frozen(DESIGN_SEAL,DESIGN_SEAL_SHA)
    artifacts=sealed['review_artifact_sha256']
    if len(artifacts)!=6:raise RuntimeError('EXACT_SIX_SOURCE_DESIGN_ARTIFACTS_REQUIRED')
    for path,digest in artifacts.items():
        if regular_sha(path)!=digest:raise RuntimeError('SOURCE_DESIGN_SEAL_ARTIFACT_CHANGED')
    load_frozen(LEDGER,LEDGER_SHA)
    candidate=design['separate_authenticated_candidate'];q=design['reference_only_QC'];contract=dict(design['resource_contract'])
    original_root=Path(contract['source_namespace']);contract['source_namespace']=str(original_root.with_name('canonical_signed_reference_source_v2'))
    contract['source_epoch_namespaces']=[str(original_root),contract['source_namespace']]
    if candidate['bytes']!=288277344 or candidate['upstream_md5']!='md5:a7773ab485827b533cb300c76356d76b' or design['archive_members']['regular_file_count']!=66 or contract['source_aggregate_cap_bytes']!=3*(1<<30):raise RuntimeError('FROZEN_SOURCE_SELECTION_OR_BUDGET_CHANGED')
    profile=P/'source_provenance/core_declared_harmonization_runtime_candidate_v4.json';profile_sha='26014a05a9ad60d42e250e6a61ad47cfb3d01eacad0db84c3d888cc88e919b57'
    A=load_frozen(profile,profile_sha)
    historical=P/'manifests/native_canonical_calibration_plan_v4_8.json';histsha='c9ad3842526a61ab24ff8fd39222c408fb611ba2bbf26e7ef0959f1e86d8de51';h=load_frozen(historical,histsha)
    base=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas');reader=base/'ldsc/ldscore/ldscore.py'
    source_pin='df98f324ca775ece05a6e8f93e0584c88bde9535c360fda6ff96bb7f475e726b'
    if regular_sha(reader)!=source_pin:raise RuntimeError('PINNED_MAINTAINED_DECODER_SOURCE_CHANGED')
    B=base/'.ldsc-env';Bsha=h['dependencies_sha256'][str(B/'bin/python')]
    if Bsha!='7dffb088cd3027e48f0127ced6f206e06111abed3b9457e580b6ebbf593c10ba':raise RuntimeError('SEALED_ORIGINAL_B_INTERPRETER_IDENTITY_CHANGED')
    Bprofile=P/'source_provenance/signed_reference_decoderB_runtime_identity_v1.json';Bprofile_sha='e1e9755f90c306aad9294b6066168ceab4e366b30038df99fb24b963d5bfea6c'
    Bmeta=load_frozen(Bprofile,Bprofile_sha)
    if Bmeta['environment_root']!=str(B) or Bmeta.get('regular_file_count')!=26334 or Bmeta.get('symlink_directory_count')!=2 or Bmeta.get('all_before_after_identities_unchanged') is not True:raise RuntimeError('ROOT_PRODUCED_B_PROFILE_METADATA_DIFFERS')
    panel=P/'statistical_validation/canonical_ld_feasibility_sources_v4/integrated_call_samples_v3.20130502.ALL.panel'
    references=dict(q['reference_sha256']);references[q['coordinate_map']['path']]=q['coordinate_map']['sha256'];references[str(panel)]='b4023dc6ee2d62ee89c8d4d347db4d348e65518d66d346574cdae7a4bbd76858'
    files=['signed_reference_source_prepare_v2.py','signed_reference_source_common_v2.py','signed_reference_source_controller_v2.py','signed_reference_source_worker_v2.py','signed_reference_decode_v2.py','signed_reference_ldsc_decoderB_v2.py','signed_reference_compact_v2.py','canonical_calibration_common_v4_8.py','30_prepare_and_run_ssd_native_campaign.py','terminal_commit_common_v2.py','extension_replay_common_v3.py']
    deps={str(P/'scripts'/name):regular_sha(P/'scripts'/name) for name in files}
    deps.update(artifacts);deps.update({str(DESIGN_SEAL):DESIGN_SEAL_SHA,str(LEDGER):LEDGER_SHA,str(profile):profile_sha,str(historical):histsha,str(reader):source_pin,str(Bprofile):Bprofile_sha})
    for path in [Path(candidate['metadata_path']),Path(design['preserved_full_protocol']['path']),P/'logs/closed_sensitivity_storage_reservation_release_v4_5.json']:
        deps[str(path)]=regular_sha(path)
    release=P/'logs/closed_sensitivity_storage_reservation_release_v4_5.json'
    closed=load_frozen(release,'ec0362b6edb45d7cb7d4cde701d9e6a1c27c6eaa684d3743171e57a9e4b28566')
    for path,digest in closed['completed_stage_evidence_sha256'].items():
        if regular_sha(path)!=digest:raise RuntimeError('EXACT_CLOSED_SENSITIVITY_COMPLETION_METADATA_CHANGED')
        deps[path]=digest
    deps['/usr/bin/curl']='5ab042572ea0e068644e3b8f9e8dd1ad197bfcf33d199316615b46ddc4390a41'
    # Existing binary and reference identities are authenticated through the
    # sealed design/old exact plan, inherited only; no body read at preparation.
    deps[str(B/'bin/python3.9')]=Bsha;deps.update(references)
    deferred={**references,str(B/'bin/python3.9'):Bsha,'/usr/bin/curl':deps['/usr/bin/curl']}
    for path,digest in deps.items():
        if path not in deferred and regular_sha(path)!=digest:raise RuntimeError('SOURCE_PREPARATION_CODE_METADATA_CHANGED')
    author=P/'reviews/signed_reference_source_author_preparation_seal_v1.json';author_sha='798be8797263ab9d8d4f4f05cf2c0966ed0d93c09dd40508fca0878f82c47311'
    old_author=load_frozen(author,author_sha)
    reject=P/'reviews/independent_signed_reference_source_prelaunch_review_v1_seal.json';reject_sha='4b024a5f1583b5d7db08f1b20b0cad6f211946a83dde9788d93af6f501b78e8b'
    old_reject=load_frozen(reject,reject_sha)
    if old_reject.get('verdict')!='REJECT_V1_OPERATIONAL_PRELAUNCH_PRESERVE_NO_EXECUTION_ADMISSION':raise RuntimeError('EXACT_V1_OPERATIONAL_REJECTION_REQUIRED')
    for closure in [old_author['file_sha256'],old_reject['file_sha256']]:
        for path,digest in closure.items():
            if path in deferred or regular_sha(path)!=digest:raise RuntimeError('OLD_V1_AUTHOR_REJECTION_CLOSURE_CHANGED')
            if path in deps and deps[path]!=digest:raise RuntimeError('OLD_NEW_DEPENDENCY_CONTRADICTION')
            deps[path]=digest
    for path,digest in {str(author):author_sha,str(reject):reject_sha,
        str(P/'logs/signed_reference_source_author_seal_root_consumed_v1.json'):'dd3f668873987c0512740eb2118dc34dbda845cb682da87d0442289a5bcd978e',
        str(P/'logs/signed_reference_source_v1_prelaunch_rejection_root_consumed_v1.json'):'eb360b19053c344afa7d6f416a7f5336213bf6f2624a0c65b4dd83b6d28f6584'}.items():
        if regular_sha(path)!=digest:raise RuntimeError('OLD_V1_ROOT_CONSUMPTION_CHANGED')
        deps[path]=digest
    numeric=sorted(q['reference_sha256'],key=lambda path:int(Path(path).name.split('.')[0]))
    if len(numeric)!=22:raise RuntimeError('EXACT22_REFERENCE_IDENTITIES_REQUIRED')
    root=Path(contract['source_namespace']);rawlock=root.parent/'extension_raw_acquisition_family.lock'
    plan=dict(schema='signed_reference_source_execution_plan_v2',prepared_utc=datetime.now(timezone.utc).isoformat(),
      scope='SOURCE_AUTHENTICATION_REFERENCE_QC_OPERATOR_ONLY',preparation_only=True,execution_admitted=False,
      output_root=str(root),self_path=str(PLAN),source_worker=str(P/'scripts/signed_reference_source_worker_v2.py'),
      source_controller=str(P/'scripts/signed_reference_source_controller_v2.py'),source_python=A['versions']['executable'],
      runtime_A_environment_root=A['environment_root'],runtime_A_profile_path=str(profile),runtime_A_profile_sha256=profile_sha,
      runtime_B_environment_root=str(B),runtime_B_required_profile_path=str(Bprofile),runtime_B_profile_sha256=Bprofile_sha,
      runtime_B_complete_current_identity='ROOT_PRODUCED_CURRENT_COMPLETE_PROFILE_BOUND;AUTHOR_READ_METADATA_ONLY;SEPARATE_INDEPENDENT_PROFILE_REVIEW_AND_FUTURE_CURRENT_FULL_HASH_GATE_REQUIRED_BEFORE_SOURCE_OPERATIONS',
      runtime_B_required_versions=dict(python='3.9.23',numpy='1.21.5',pandas='1.3.3',bitarray='2.8.3'),
      decoderB_python=str(B/'bin/python'),decoderB_wrapper=str(P/'scripts/signed_reference_ldsc_decoderB_v2.py'),
      maintained_decoderB_path=str(reader),maintained_decoderB_sha256=source_pin,
      decoderB_contract='UNCHANGED_ORIGINAL_PLINKBEDFILE_NEXTSNPS;<=288_PRESPECIFIED_ALREADY_ELIGIBLE_SNP_ROWS;ASSERT_ORIGINAL_FILTER_REMOVES_NONE;A1=2-A2;TARGET_G=-Y_B*FROZEN_SIGN;ORIGINAL_PADDING_AND_ZERO_FALLBACK_NO_INDEPENDENT_QC_CREDIT',
      compact_representation=SCHEMA,compact_schema_sha256=SCHEMA_SHA256,
      conservative_compact_capacity_qualification='COMPLETE_RETENTION_WITHIN_FIXED_STRUCTURAL_BUDGET;SOURCE_DEPENDENT_UTF8_POOL_AND_SQLITE_LIMITS_FAIL_WHOLE_STAGE_BEFORE_GENOTYPE_QC;NO_SOURCE_FIT_CERTIFICATION',
      prior_v1_author_seal_sha256=author_sha,prior_v1_rejection_seal_sha256=reject_sha,
      fixed_directory_contract=['','archive','extracted','extracted/members','qc','qc/controls','tmp'],
      curl_contract='FIRST_POSITION_-q_DISABLE_DEFAULT_CURLRC;EXACT_HTTPS200_NO_REDIRECT_RESUME_RETRY',
      candidate=candidate,archive_members=design['archive_members'],reference_only_QC=q,independent_controls=design['independent_controls'],
      finite_reference_truth_contract=design['finite_reference_truth_contract'],preserved_full_protocol=design['preserved_full_protocol'],
      expected_master_reference_rows=1290028,reference_files_numeric_order=numeric,coordinate_map_path=q['coordinate_map']['path'],
      authoritative_panel_path=str(panel),reference_body_sha256=references,dependency_sha256=deps,
      preparation_deferred_dependency_sha256=deferred,preparation_deferred_qualification='SEALED_INHERITED_IDENTITIES_NOT_CURRENT_BODY_HASHES;MANDATORY_FULL_CURRENT_GATE_AT_FUTURE_INDEPENDENTLY_ADMITTED_EXECUTION',
      resource_contract=contract,closed_sensitivity_namespace=str(root.parent/'sensitivities'),
      raw_transfer_family_mutex=str(rawlock),owned_monitor_helper=str(P/'scripts/30_prepare_and_run_ssd_native_campaign.py'),curl='/usr/bin/curl',
      worker_count=1,worker_command_template=[A['versions']['executable'],'-B',str(P/'scripts/signed_reference_source_worker_v2.py'),'--plan',str(PLAN),'--plan-sha256','FIXED_REVIEWED_PLAN_SHA','--started-monotonic','CONTROLLER_STAGE_START',
        '--admission',str(ADMISSION),'--admission-sha256','FIXED_REVIEWED_ROOT_ADMISSION_SHA','--heavy-lock-fd','INHERITED_SHARED_HEAVY_FD','--raw-lock-fd','INHERITED_RAW_FAMILY_FD'],
      G1='OPTIONAL_OLD_LOCAL_CONTINUITY_NOT_REQUESTED_NOT_RUN_NO_OLD_COMPACT_AUTHENTICATION_CREDIT',
      G2_G3='NEW_AUTHENTICATED_CANDIDATE_EXACT_MD5_BYTES_MEMBER_SAFETY_FAM_BIM_BED;NO_RETRY_RESUME_MIRROR',
      G4='FREEZE_ALL_HISTORICAL_MASTER_RANKS_STATIC_EXACT_GRCh37_COORDINATE_AND_ALLELE_INTERSECTION_BEFORE_GENOTYPE_QC;NO_GWAS',
      G5='POST_J_PRESPECIFIED_WINDOWS_SCATTER_MAINTAINED_B_INDEPENDENCE_SIGNED_GRAM_AND_FULL_CHROMOSOME_STREAMING_OPERATORS',
      source_build_qualification='PUBLISHER_ARCHIVE_NAME_DOES_NOT_ESTABLISH_FULL_SOURCE_BUILD;ONLY_EXACT_STATIC_SELECTED_GRCh37_COMPATIBILITY_ALLOWED;NO_UNRESOLVED_SELECTED_COORDINATE_ENTERS_QC',
      private_data_policy='ALL_ACTUAL_FAM_IDS_SNP_GENOTYPE_QC_AND_INDEPENDENT_EXPORTS_REMAIN_PRIVATE_SSD_NAMESPACE;REPOSITORY_AGGREGATE_HASH_RECEIPTS_ONLY',
      terminal_contract='PRIVATE_IMMUTABLE_OWNED_TERMINAL2_PENDING;INTENDED_PAYLOADS_BEFORE_WRITES;IMMUTABLE_BEFORE_AFTER_PARSE;EXACT_REGISTERED_APPLEDOUBLE;FULL_POSTPERSIST_IDENTITY_RUNTIME_RESOURCE_SIGNAL_AND_CLEANUP_GATES',
      no_source_operation_until_separate_independent_review_and_root_admission=True,
      S_simulations_run=0,C_calibration_credit=0,estimator_calls=0,biological41Cov_release_admitted=False,
      remaining_launch_prerequisites=['Separate independent complete source/controller/decoder code and actual-plan review','Separate independent review of the bound current B readonly runtime profile','Root exact source-only admission and current physical SSD/resource/dual-lock preflight','No active100/13 acquisition or heavy worker; both exact shared locks must acquire nonblocking','Source packaging/build/sample/QC/window/operator gates succeed on actual authenticated bytes; no pre-execution achievement claim'])
    digest=save(PLAN,plan)
    print(json.dumps(dict(status='PREPARATION_ONLY_NO_SOURCE_BODY_OPERATION',plan=str(PLAN),sha256=digest,
      deferred_dependency_count=len(deferred),workers_launched=0,source_downloads=0,execution_admitted=False)))


if __name__=='__main__':prepare()

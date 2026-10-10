from pathlib import Path
import sys, json, importlib.util, time, shutil, resource, os, hashlib
from datetime import datetime, timezone

P=Path('/Users/swethasunilkumar/Documents/Codex/2026-10-08/go/work/sleep-gwas-atlas/sleep_unified_research_v4')
sys.path.insert(0,str(P/'scripts'))
from extension_replay_common_v5 import sha, body_hashes, acquisition_family_gate, acquisition_operational_gate, acquisition_receipt_gate, physical_mount
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
PLAN=SSD/'extension_pipeline_replay_v7/extension_pipeline_replay_plan_v7.json'
PLAN_SHA='ad0c08bc0b0a62c8a1bc070221689f2733b0eefe8f25de8cc1bd8689a81b1ba5'
ADMISSION=P/'manifests/extension_pipeline_replay_admission_v7.json'
SEAL=P/'reviews/independent_extension_pipeline_v7_prelaunch_review_seal_v1.json'
SEAL_SHA='56a857ff7440382c133ca97362a7046486d8924ceccf7ee4d2737df5622622c8'
PROOF=P/'logs/extension_pipeline_v7_root_consumed_independent_prelaunch_review_v1.json'
assert not ADMISSION.exists() and not ADMISSION.is_symlink()
assert sha(PLAN)==PLAN_SHA and sha(SEAL)==SEAL_SHA
plan=json.loads(PLAN.read_text());seal=json.loads(SEAL.read_text())
assert len(plan['dependencies_sha256'])==518 and len(plan['archived_input_sha256'])==100
assert seal['verdict']=='QUALIFIED_PREFLIGHT_PASS'
assert sha(PROOF)=='61f19a29df19ead9898858fe289d0be4d7724221746f0610d9e0a9f8a4ce6178'
proof=json.loads(PROOF.read_text())
assert proof['all_files_and_seal_current'] and proof['verified_regular_metadata_file_count']==39522
review_map={name:digest for name,digest in seal['file_sha256'].items()
    if Path(name).parent==P/'reviews' and Path(name).name.startswith('independent_extension_pipeline_v7_')}
review_map.update({str(SEAL):SEAL_SHA,str(PROOF):sha(PROOF)})
assert any(name.endswith('.md') for name in review_map) and any(name.endswith('.json') for name in review_map)
for name,digest in review_map.items():assert sha(name)==digest,name
started=time.monotonic()

def state():
    internal=shutil.disk_usage('/System/Volumes/Data').free
    free=shutil.disk_usage(SSD).free
    assert internal>=plan['guard']['internal_floor_bytes'] and free>=plan['guard']['SSD_floor_bytes']
    assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<256<<20 and time.monotonic()-started<7200
    return dict(internal_free_bytes=internal,SSD_free_bytes=free)

state();mount=physical_mount();acquisition_operational_gate(plan)
family=acquisition_family_gate(plan,require_complete=True)
assert len(family['source_receipt_sha256'])==100
profile_path=P/'source_provenance/signed_reference_decoderB_runtime_identity_v1.json'
profile_sha='e1e9755f90c306aad9294b6066168ceab4e366b30038df99fb24b963d5bfea6c'
profile_seal=P/'reviews/independent_decoderB_runtime_identity_review_v1_seal.json'
profile_proof=P/'logs/decoderB_runtime_independent_review_root_consumed_v1.json'
assert sha(profile_path)==profile_sha
assert sha(profile_seal)=='7d7f388d4380632682f3dde920b32f93e2b7e454d90b2a72dac1afde3c7a0a51'
assert sha(profile_proof)=='10f53c18adc86bb78b743b43951460b51050e832a5fc449812517bc897c16441'
profile=json.loads(profile_path.read_text())
logical_python='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas/.ldsc-env/bin/python'
logical_reference='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas/work/track_b_completion/local_dependency_copies/ref/w_hm3.snplist'
assert plan['w_hm3']==logical_reference
allowed_file_links={
 logical_python:dict(literal_link='python3.9',resolved_target=str(Path(logical_python).parent/'python3.9'),target_is_file=True,resolved_sha256='7dffb088cd3027e48f0127ced6f206e06111abed3b9457e580b6ebbf593c10ba'),
 logical_reference:dict(literal_link='eur_w_ld_chr/w_hm3.snplist',resolved_target=str(Path(logical_reference).parent/'eur_w_ld_chr/w_hm3.snplist'),target_is_file=True,resolved_sha256='ec73fca0b696e8beba465b51e52911676fcade375bcc9475d99c0ec30509d3ed')}
assert profile['symlinks'][logical_python]==allowed_file_links[logical_python]
current_deps={};dependency_bytes=0;current_declared_file_links={}
for name,digest in plan['dependencies_sha256'].items():
    q=Path(name)
    assert q.is_file() and not any(parent.is_symlink() for parent in q.parents),name
    if q.is_symlink():
        assert name in allowed_file_links,name
        route=allowed_file_links[name]
        assert os.readlink(q)==route['literal_link'] and str(q.resolve(strict=True))==route['resolved_target'],name
        target=Path(route['resolved_target'])
        assert target.is_file() and not target.is_symlink() and not any(parent.is_symlink() for parent in target.parents),name
        assert digest==route['resolved_sha256'] and sha(target)==digest,name
        current_declared_file_links[name]=dict(route,bytes=target.stat().st_size)
    else:
        assert name not in allowed_file_links,name
    assert sha(q)==digest,name
    current_deps[name]=digest;dependency_bytes+=q.stat().st_size
assert set(current_declared_file_links)==set(allowed_file_links)

state()
print(json.dumps(dict(verified_current_dependencies=518,archived_inputs=100,dependency_bytes=dependency_bytes)),flush=True)
body_map={};source_map={}
for index,member in enumerate(plan['members'],1):
    receipt=Path(member['acquisition_receipt']);digest=family['source_receipt_sha256'][str(receipt)]
    acquisition_receipt_gate(plan,member,digest);source_map[str(receipt)]=digest
    q=Path(member['acquisition_member']['body_path'])
    assert not any(parent.is_symlink() for parent in q.parents),str(q)
    got=body_hashes(q);expected=member['acquisition_member']
    assert got==dict(sha256=expected['expected_sha256'],md5=expected['expected_md5'],bytes=expected['expected_bytes']),str(q)
    assert sha(receipt)==digest
    body_map[str(q)]=got;state()
    if index%5==0:print(json.dumps(dict(current_full_compressed_body_identities_verified=index,of=100)),flush=True)
assert acquisition_family_gate(plan,source_map,require_complete=True)==family
spec=importlib.util.spec_from_file_location('root_extension7_pinned_owned_gate',P/'scripts/sensitivity_executor_v4_4.py')
protected=importlib.util.module_from_spec(spec);spec.loader.exec_module(protected)
historical=protected.baseline_gate(plan);state();physical_mount()
assert sha(PLAN)==PLAN_SHA and sha(SEAL)==SEAL_SHA
for name,digest in review_map.items():assert sha(name)==digest,name
namespace_bytes=sum(q.stat().st_size for q in SSD.rglob('*') if q.is_file() and not q.is_symlink())
assert namespace_bytes<=plan['guard']['global_reservation_bytes']
spec=importlib.util.spec_from_file_location('root_extension7_exact_controller',P/'scripts/47_run_extension_pipeline_replay_v7.py')
controller=importlib.util.module_from_spec(spec);spec.loader.exec_module(controller)
controller.no_preserved_outputs_gate(plan)
assert not (controller.OUT/'pipeline_attempt_v7.json').exists()
binding=plan['acquisition_execution_binding']
admission=dict(schema='root_actual_extension7_historical100_preprocessing_admission_v2',
    recorded_utc=datetime.now(timezone.utc).isoformat(),execution_admitted=True,
    plan_sha256=PLAN_SHA,executor_sha256=sha(P/'scripts/47_run_extension_pipeline_replay_v7.py'),
    checkpoint_receipt_binding_policy=plan['checkpoint_receipt_binding_policy'],
    acquisition_receipt_sha256=source_map,acquisition_execution_plan_sha256=plan['acquisition_plan_sha256'],
    acquisition_executor_sha256=binding['executor_sha256'],
    acquisition_root_admission_sha256=binding['sha256'][binding['root_admission']],
    independent_review_sha256=review_map,independent_actual_plan_review_pass=True,
    all518_current_dependency_sha256=current_deps,all100_current_compressed_body_identity=body_map,
    current_declared_file_symlink_routes=current_declared_file_links,
    current_decoderB_runtime_identity_profile_sha256=profile_sha,
    independently_reviewed_decoderB_runtime_identity_root_proof_sha256=sha(profile_proof),
    preserved_root_checker_v1_stop_receipt_sha256=sha(P/'logs/extension_pipeline_v7_root_identity_checker_v1_stop_receipt.json'),
    full100_acquisition_terminal_evidence=family,historical190_gate_receipt_sha256=historical,
    current_resource_state=state(),global_campaign_regular_non_symlink_bytes=namespace_bytes,
    dependency_bytes_verified=dependency_bytes,
    compressed_raw_body_bytes_verified=sum(value['bytes'] for value in body_map.values()),
    workers_launched_at_admission=0,estimator_calls=0,automatic_retry=False,
    scientific_membership_or_threshold_changes=False,
    current_GZIP_EOF_or_ordered_preprocessing_content_verified=False,
    independent_replication_or_source_permission_admitted=False,
    launch_schedule='Queued until admitted standalone FinnGen/canonical/source stages can use the shared heavy lock after the active core35 campaign; no concurrent scientific worker or modified lock protocol.',
    qualification='Qualified original100 raw-to-template processing admission after genuine complete full100 producer and fresh actual input/dependency/output identities. Entire ordered harmonized/munged content, literal QC and gzip CRC/EOF still require real execution and independent result review. Original filter-count concordance remains unestablished when historical QC TSV is unavailable.',
    root_identity_elapsed_seconds=time.monotonic()-started,
    peak_root_identity_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    root_current_input_identity_checker_path=str(Path(__file__)),
    root_current_input_identity_checker_sha256=sha(__file__))
assert admission['executor_sha256']=='6abd57c9d743ea1c33b7ccbed0f88c51a7af0972db24655284f3e279fa483f00'
controller.admission_gate(plan,admission,PLAN_SHA,admission['executor_sha256'])
payload=(json.dumps(admission,indent=2,allow_nan=False)+'\n').encode()
intended=hashlib.sha256(payload).hexdigest()
with ADMISSION.open('xb') as stream:stream.write(payload);stream.flush();os.fsync(stream.fileno())
assert sha(ADMISSION)==intended
print(json.dumps(dict(status='ACTUAL_EXTENSION7_ROOT_ADMITTED_QUEUED',admission_path=str(ADMISSION),
    admission_sha256=intended,current_dependencies=518,full_current_source_bodies=100,
    historical_receipts=len(historical),elapsed_seconds=admission['root_identity_elapsed_seconds'])),flush=True)

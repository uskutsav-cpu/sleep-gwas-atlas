"""Source-only identities/resources. Execution is separately admitted; no fits."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import time
from contextlib import contextmanager
from canonical_calibration_common_v4_8 import (regular_sha, intended_sha, write_intended_new,
    terminal_intended_sha, freeze_companions, deferred_termination_signals,
    assert_no_termination, exclusive_heavy_lock, safe_diagnostic, TERMINATION_REQUEST)

P=Path(__file__).resolve().parents[1]
PLAN=P/'manifests/signed_reference_source_execution_plan_v1.json'
ADMISSION=P/'manifests/signed_reference_source_execution_admission_v1.json'
DESIGN=P/'manifests/signed_ld_503eur_source_only_prelaunch_plan_v4_1.json'
DESIGN_SHA='54d62c6361ae57fb8791e1b55271afc1b806088902c52ddfbca524b7511847be'
DESIGN_SEAL=P/'reviews/signed_ld_503eur_source_design_seal_v4_1.json'
DESIGN_SEAL_SHA='e63162e196fecece21f6114fe22d2563ff0857f61efb5ec4972af9c7cef60ead'
LEDGER=P/'manifests/global_SSD_resource_reservation_v4_5.json'
LEDGER_SHA='7ce69095a6f67ae681d547476d1557afac2de0f97d083b12bc1b689d0788417f'
INTENDED_OUTPUT_SHA256={}


def record_output(path,expected):
    key=str(Path(path))
    if key in INTENDED_OUTPUT_SHA256 and INTENDED_OUTPUT_SHA256[key]!=expected:raise RuntimeError('SOURCE_FIRST_OUTPUT_IDENTITY_REBOUND')
    if regular_sha(path)!=expected:raise RuntimeError('SOURCE_OUTPUT_EXPECTED_HASH_MISMATCH')
    INTENDED_OUTPUT_SHA256.setdefault(key,expected)
    return expected


def load_frozen(path, expected):
    if regular_sha(path)!=expected:raise RuntimeError('SOURCE_METADATA_CHANGED_BEFORE_PARSE: '+str(path))
    value=json.loads(Path(path).read_text())
    if regular_sha(path)!=expected:raise RuntimeError('SOURCE_METADATA_CHANGED_AFTER_PARSE: '+str(path))
    return value


def save(path, value):
    expected=intended_sha(value);write_intended_new(path,value,expected)
    return record_output(path,expected)


def binary_new(path,payload):
    path=Path(path);expected=hashlib.sha256(payload).hexdigest()
    if any(p.is_symlink() for p in path.parents):raise RuntimeError('SOURCE_OUTPUT_PARENT_SYMLINK')
    with path.open('xb') as f:f.write(payload);f.flush();os.fsync(f.fileno())
    fd=os.open(path.parent,os.O_RDONLY)
    try:os.fsync(fd)
    finally:os.close(fd)
    if regular_sha(path)!=expected:raise RuntimeError('SOURCE_INTENDED_BINARY_CHANGED')
    return record_output(path,expected)


def streamed_new(path, records, guard=lambda:None):
    """Freeze intended bytes as they are generated, before any first read."""
    path=Path(path)
    if any(p.is_symlink() for p in path.parents):raise RuntimeError('SOURCE_OUTPUT_PARENT_SYMLINK')
    h=hashlib.sha256();count=0
    with path.open('xb') as f:
        for payload in records:
            guard();h.update(payload);f.write(payload);count+=1
        f.flush();os.fsync(f.fileno())
    fd=os.open(path.parent,os.O_RDONLY)
    try:os.fsync(fd)
    finally:os.close(fd)
    expected=h.hexdigest()
    if regular_sha(path)!=expected:raise RuntimeError('SOURCE_STREAMED_INTENDED_BYTES_CHANGED')
    record_output(path,expected)
    return expected,count


def full_tree_bytes(root):
    root=Path(root)
    if not root.exists():return 0
    if root.is_symlink() or any(p.is_symlink() for p in root.parents):raise RuntimeError('SOURCE_METER_ROOT_SYMLINK')
    total=0
    for path in root.rglob('*'):
        if path.is_symlink():raise RuntimeError('SOURCE_CURRENT_METER_SYMLINK: '+str(path))
        if path.is_file():total+=path.stat().st_size
        elif not path.is_dir():raise RuntimeError('SOURCE_CURRENT_METER_SPECIAL_FILE')
    return total


def resource_snapshot(plan,started,launch=False):
    root=Path(plan['output_root']);contract=plan['resource_contract']
    total=full_tree_bytes(root);a=full_tree_bytes(root/'archive');e=full_tree_bytes(root/'extracted')
    q=total-a-e
    outside=root.with_name('._'+root.name)
    if outside.exists() or outside.is_symlink():
        regular_sha(outside);q+=outside.stat().st_size;total+=outside.stat().st_size
    campaign=full_tree_bytes(root.parent);closed=full_tree_bytes(plan['closed_sensitivity_namespace'])
    state={'internal_free_bytes':shutil.disk_usage('/System/Volumes/Data').free,
      'ssd_free_bytes':shutil.disk_usage(root.parent).free,'source_aggregate_bytes':total,
      'archive_bytes':a,'extraction_bytes':e,'QC_controls_metadata_bytes':q,
      'global_campaign_bytes':campaign,'closed_sensitivity_bytes':closed,'elapsed_seconds':time.monotonic()-started}
    failures=[]
    if state['internal_free_bytes']<contract['internal_floor_bytes']:failures.append('INTERNAL3GIB_FLOOR')
    floor=contract['SSD_launch_free_required_bytes'] if launch else contract['SSD_floor_bytes']
    if state['ssd_free_bytes']<floor:failures.append('SSD8GIB_LAUNCH_OR5GIB_RUNTIME_FLOOR')
    if total>contract['source_aggregate_cap_bytes'] or a>contract['archive_and_all_failed_partial_cap_bytes'] or e>contract['extracted_all_partial_and_failed_cap_bytes'] or q>contract['QC_controls_metadata_and_all_failed_output_cap_bytes']:failures.append('SOURCE3GIB_OR_SUBCAP')
    if campaign>contract['global_ceiling_bytes']:failures.append('GLOBAL300GIB_CURRENT_METER')
    if closed>contract['closed_sensitivity_cap_bytes']:failures.append('CLOSED_SENSITIVITY256MIB_RETENTION')
    if state['elapsed_seconds']>contract['whole_stage_monotonic_deadline_seconds']:failures.append('WHOLE_STAGE2H_DEADLINE')
    assert_no_termination()
    if failures:raise RuntimeError('SOURCE_RESOURCE_GUARD: '+','.join(failures))
    return state


class Guard:
    def __init__(self,plan,started):self.plan=plan;self.started=started;self.last=-1
    def __call__(self,force=False):
        assert_no_termination()
        if time.monotonic()-self.started>7200:raise RuntimeError('WHOLE_SOURCE_STAGE_DEADLINE')
        if force or time.monotonic()-self.last>=1:
            resource_snapshot(self.plan,self.started);self.last=time.monotonic()


def check_dependencies(plan,include_reference=True):
    deferred=plan['reference_body_sha256'] if not include_reference else {}
    for path,digest in plan['dependency_sha256'].items():
        if path not in deferred and regular_sha(path)!=digest:raise RuntimeError('SOURCE_DEPENDENCY_CHANGED: '+path)


def identity_gate(plan,plan_path,admission_path,plan_sha,admission_sha):
    if regular_sha(plan_path)!=plan_sha or regular_sha(admission_path)!=admission_sha:raise RuntimeError('SOURCE_FIXED_PLAN_OR_ADMISSION_CHANGED')
    a=load_frozen(admission_path,admission_sha)
    if a.get('execution_admitted') is not True or a.get('scope')!='SOURCE_AUTHENTICATION_REFERENCE_QC_OPERATOR_ONLY' or a.get('plan_sha256')!=plan_sha or a.get('independent_source_code_review_pass') is not True or a.get('independent_decoderB_contract_pass') is not True:raise RuntimeError('NO_INDEPENDENT_SOURCE_EXECUTION_ADMISSION')
    if a.get('allow_S_simulation') is not False or a.get('allow_estimators') is not False or a.get('allow_C_calibration_credit') is not False or a.get('allow_biological41Cov') is not False:raise RuntimeError('SOURCE_ONLY_SCOPE_REQUIRED')
    if not a.get('review_artifact_sha256'):raise RuntimeError('SOURCE_REVIEW_MAP_REQUIRED')
    for path,digest in a['review_artifact_sha256'].items():
        if regular_sha(path)!=digest:raise RuntimeError('SOURCE_INDEPENDENT_REVIEW_CHANGED')
    check_dependencies(plan)
    if regular_sha(LEDGER)!=LEDGER_SHA:raise RuntimeError('SOURCE_LEDGER_CHANGED')
    from extension_replay_common_v3 import physical_mount
    physical_mount()


@contextmanager
def shared_raw_family_lock(path,before_release):
    path=Path(path)
    if not path.is_file() or path.is_symlink() or any(p.is_symlink() for p in path.parents):raise RuntimeError('EXACT_EXISTING_RAW_FAMILY_LOCK_REQUIRED')
    with path.open('r+b') as f:
        fcntl.flock(f.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:yield f.fileno()
        finally:
            while True:
                try:before_release();break
                except BaseException as error:
                    safe_diagnostic(lambda:'SOURCE_RAW_FAMILY_LOCK_RETAINED_UNVERIFIED_CLEANUP: '+repr(error))
                    try:time.sleep(1)
                    except BaseException:pass
            fcntl.flock(f.fileno(),fcntl.LOCK_UN)


def verify_outputs(mapping):
    for path,digest in mapping.items():
        if regular_sha(path)!=digest:raise RuntimeError('FROZEN_SOURCE_OUTPUT_CHANGED: '+path)


def freeze_registered(root,targets,frozen):
    """Strict registered ordinary outputs plus exact regular AppleDouble only."""
    targets=[Path(p) for p in targets]
    ordinary={str(p) for p in targets if p.is_file() and not p.name.startswith('._')}
    for p in ordinary:
        digest=regular_sha(p)
        if p in frozen and frozen[p]!=digest:raise RuntimeError('SOURCE_REGISTERED_OUTPUT_REBOUND')
        frozen.setdefault(p,digest)
    companion_map={p:h for p,h in frozen.items() if Path(p).name.startswith('._')}
    freeze_companions(targets,companion_map)
    frozen.update(companion_map);verify_outputs(frozen)
    observed=set()
    for path in Path(root).rglob('*'):
        if path.is_symlink() or (not path.is_dir() and not path.is_file()):raise RuntimeError('SOURCE_OUTPUT_SYMLINK_OR_SPECIAL')
        if path.is_file():observed.add(str(path))
    outside=Path(root).with_name('._'+Path(root).name)
    if outside.exists() or outside.is_symlink():observed.add(str(outside))
    if observed!=set(frozen):raise RuntimeError('SOURCE_EXACT_REGISTERED_OUTPUT_INVENTORY_DIFFERS')

"""Targeted stopped-v3 metadata and exact90 meter-function witness only.

No source body or successful11 receipt reread, real campaign traversal, download,
worker, process probe or old suite. The preserved race witness uses a private
one-byte file removed after is_file by a local wrapper, not a running core spool.
"""
import ast
import hashlib
import json
from pathlib import Path
import resource
import time

P=Path(__file__).resolve().parents[1]
R=P/'reviews'
S=P/'scripts'
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
V3=SSD/'validation_raw_replay_v3'
FAMILY=V3/'validation_raw_acquisition_family_receipt_v4_3.json'
SOURCE=V3/'receipts/finngen_r13_K11_CHOLELITH.json'
CODE=S/'90_acquire_validation_raw_sources_v3.py'
MONITOR=S/'30_prepare_and_run_ssd_native_campaign.py'
OUT=R/'independent_validation13_resource_race_receipt_v1.json'
FIXTURES=R/'independent_validation13_resource_race_fixtures_v1'
started=time.monotonic()
fixed={
    str(FAMILY):'995a00139dfbb40f772dfdd262ad93f20fc9fbef65acddf9b2cb50d46e3a6a3d',
    str(P/'logs/validation_raw_acquisition_family_receipt_v4_3.json'):'995a00139dfbb40f772dfdd262ad93f20fc9fbef65acddf9b2cb50d46e3a6a3d',
    str(SOURCE):'bde5ece0dfb1aa799502a7f57ca57b6856e7c0130fd7ccbf3964784ec73ae662',
    str(P/'source_provenance/validation_raw_acquisition_v4_3/finngen_r13_K11_CHOLELITH.json'):'bde5ece0dfb1aa799502a7f57ca57b6856e7c0130fd7ccbf3964784ec73ae662',
    str(V3/'family_pending_v4_3.json'):'afe2a4855931c4c7f42641c04ba609a5d41aeeb067457b81573ebff43763f6fb',
    str(V3/'family_ownership.json'):'01e0df53ff81c7d2d6e8e16e9bbb847b1a3ab8ac70260d08fe39180b2f86ba92',
    str(CODE):'8d0765930a346d1097d15a7b5fc3bea294f782b4eb8f52a98f894bbf693152cc',
}


def sha(path):
    path=Path(path)
    assert path.is_file() and not path.is_symlink()
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path,value):
    with path.open('x') as stream:json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')


assert all(sha(q)==v for q,v in fixed.items())
fixed[str(MONITOR)]=sha(MONITOR)
family=json.loads(FAMILY.read_text())
source=json.loads(SOURCE.read_text())
pending=json.loads((V3/'family_pending_v4_3.json').read_text())
assert family['status']=='FAMILY_STOP_PRESERVED_REQUIRES_EXPLICIT_AUDIT'
assert len(family['source_receipts'])==family['completed_source_count']==11
assert len({x['path'] for x in family['source_receipts']})==11
assert family['plan_sha256']==source['plan_sha256']==pending['binding']['plan_sha256']=='8e93d548163bda284b6292b9780962060bbc3ff4afde44e5a916c80f3b69b708'
assert source['status']=='FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
assert source['member']['index']==12 and source['member']['source_id']=='finngen_r13_K11_CHOLELITH'
assert source['retained_partial_bytes']==57671680 and source['retained_partial_sha256']=='3537a1b14c7e3560c43d38898dfadd6ba280738498ef5de4c21195c962a1a415'
assert source['member']['expected_bytes']==809815233
assert source['member']['expected_sha256']=='50af445a7edcd53129da6e4f21a461af62ae9dd53ad983f00525bb2f29a298e1'
assert source['member']['expected_md5']=='29db68b6a610c94cb5b8d48b7b01ccb7'
assert source['resume_offset']==0 and source['returncode']==-15
assert source['teardown']['teardown_verified'] is True and source['teardown']['remaining_group_members']==[]
assert source['teardown']['signals']==['SIGTERM'] and not source['teardown'].get('cleanup_error')
assert len(source['samples'])==1
assert source['stop_reason'].startswith('FileNotFoundError:') and 'ephemeral_bounded_spool/columns/typed_08_00450000.pkl' in source['stop_reason']
assert source['stop_reason'] in family['stop_reason']
assert all(r['internal_free_bytes']>=3<<30 and r['ssd_free_bytes']>=5<<30
           for r in [source['resource_before'],source['resource_after'],family['resource_after']]+source['samples'])
assert source['elapsed_seconds']<7200 and family['elapsed_seconds']<96*3600
assert source['peak_observed_owned_rss_bytes']<2<<30
assert not (V3/'family_terminal_seal_v4_3.json').exists()
tree=ast.parse(CODE.read_text())
meter=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='global_namespace_gate')
assert "p.stat().st_size" in ast.unparse(meter) and "p.is_file()" in ast.unparse(meter)
monitor_tree=ast.parse(MONITOR.read_text())
cleanup=next(n for n in monitor_tree.body if isinstance(n,ast.FunctionDef) and n.name=='terminate_owned')
assert 'proc.wait(timeout=10)' in ast.unparse(cleanup)
FIXTURES.mkdir(exist_ok=False)
results=[]
for name,kind,expected in [('healthy','healthy','PASS'),('removed_after_is_file','removed','FileNotFoundError'),
                           ('permission_failure','permission','PermissionError'),('over_cap','healthy','AssertionError'),
                           ('symlink_excluded','symlink','PASS')]:
    folder=FIXTURES/name;folder.mkdir()
    q=folder/'tiny.txt';q.write_bytes(b'x')
    events=[]
    class Entry:
        def is_file(self):
            events.append('is_file')
            present=q.is_file()
            if kind=='removed':q.unlink();events.append('own_tiny_file_removed')
            return present
        def is_symlink(self):
            events.append('is_symlink');return kind=='symlink'
        def stat(self):
            events.append('stat')
            if kind=='permission':raise PermissionError('PRIVATE_CONTROL_PERMISSION_DENIED')
            return q.stat()
    class Walk:
        def rglob(self,pattern):
            assert pattern=='*';return [Entry()]
    ns={'SSD':Walk()}
    exec(compile(ast.Module(body=[meter],type_ignores=[]),str(CODE),'exec'),ns)
    value=None;error=None
    try:
        value=ns['global_namespace_gate']({'SSD_reservation_bytes':0 if name=='over_cap' else 300<<30})
        actual='PASS'
    except Exception as exc:
        actual=type(exc).__name__;error=str(exc)
    assert actual==expected
    if name=='healthy':assert value==1
    if name=='symlink_excluded':assert value==0 and 'stat' not in events
    case=dict(control=name,actual=actual,expected=expected,events=events,observed_bytes=value,error=error,
        uses_only_own_tiny_metadata=True,actual_campaign_meter_traversal=False)
    save(folder/'witness.json',case);results.append(case)
assert all(sha(q)==v for q,v in fixed.items())
files={str(q):sha(q) for q in FIXTURES.rglob('*') if q.is_file() and not q.is_symlink()}
assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<64<<20 and time.monotonic()-started<60
result=dict(schema='independent_validation13_v3_targeted_resource_race_v1',
    exact_failed_source_family_and_resource_code_metadata_sha256=fixed,
    recorded_v3_plan_sha256=family['plan_sha256'],recorded_admission_sha256=pending['binding']['admission_sha256'],
    completed11_source_receipt_bindings_from_failed_family_only=family['source_receipts'],
    completed11_receipts_or_body_files_reread=False,
    v3_failed_family_and_PENDING_preserved=True,source12_receipt_status=source['status'],
    source12_returncode=source['returncode'],source12_teardown=source['teardown'],
    source12_prefix_bytes=source['retained_partial_bytes'],source12_prefix_sha256_recorded_not_freshly_rehashed=source['retained_partial_sha256'],
    source12_expected_original_member=source['member'],stop_reason=source['stop_reason'],
    floor_evidence=dict(source_before=source['resource_before'],source_after=source['resource_after'],family_after=family['resource_after']),
    cause='Observed path disappearance between p.is_file and p.stat in non-atomic global namespace census; propagates ENOENT and acquisition fail-stops.',
    measured_operational_race_reproduced_with_exact_meter_AST=True,source_QC_or_identity_failure_demonstrated=False,
    actual_global_cap_at_stop_independently_measured=False,actual_core_worker_OR_live_spool_inspected=False,
    actual_resource_process_census_OR_mount_probe_performed=False,
    controls=results,control_count=5,fixture_regular_sha256=files,
    no_corrected_production_executor_or_meter_authored=True,no_resume_or_acquisition_executed=True,
    unchanged_required_caps=dict(global_bytes=300<<30,internal_floor=3<<30,ssd_floor=5<<30,
        transfer_RSS_limit=2<<30,source_seconds=7200,family_seconds=96*3600),
    source12_resume_currently_disabled=True,current_header_oracle_accepts_200_only=True,
    preferred_continuation='Distinct reviewed/admitted epoch: immutable reuse11 old authenticated sources; new full source12 and source13; preserve failedsource12 prefix and originalv3PENDING/master/receipts; frozen mixed-origin13 provenance and final full current authentication.',
    new_preserved_prefix_reservation_bytes=57671680,
    full13_original_body_declared_budget_bytes=10058648185,
    minimum_original13_plus_this_failed_prefix_bytes=10058648185+57671680,
    Range_remaining_network_body_bytes=809815233-57671680,
    Range_qualification='Only after separately frozen pinned-generation/ETag/206/exactContent-Range/currentprefixcopySHA/fulloriginalSHA_MD5_bytes+teardown/terminal oracle; no existingv3range path admission.',
    elapsed_seconds=time.monotonic()-started,max_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
save(OUT,result)
print(json.dumps(dict(diagnosis='RESOURCE_CENSUS_ENOENT_RACE_CONFIRMED',controls=5,
    family_completed=11,source12_prefix_bytes=57671680,receipt_sha256=sha(OUT),
    elapsed_seconds=result['elapsed_seconds'],max_RSS_bytes=result['max_RSS_bytes'])))

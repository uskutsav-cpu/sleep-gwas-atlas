"""Focused tiny authored delta controls; no production bodies/locks/transfers."""
import ast
import contextlib
import copy
import errno
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import resource
import sys
import time
from types import SimpleNamespace

P=Path(__file__).resolve().parents[1]
R=P/'reviews';CODE=P/'scripts/108_acquire_validation_raw_sources_v4.py'
FIX=R/'validation13_raw_continuation_author_fixtures_v4'
OUT=R/'validation13_raw_continuation_author_controls_v4.json'
started=time.monotonic()
sys.path.insert(0,str(P/'scripts'))
spec=importlib.util.spec_from_file_location('_author_validation108',CODE);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
FIX.mkdir(exist_ok=False)
results=[]
checks=0

def check(condition):
    global checks
    assert condition
    checks+=1

def case(name,call,expected=None):
    try:
        value=call()
        check(expected is None)
        results.append(dict(case=name,status='PASS',expected='PASS'))
        return value
    except Exception as exc:
        if expected is None:raise
        check(isinstance(exc,expected))
        results.append(dict(case=name,status='EXPECTED_REJECTION',error=type(exc).__name__+': '+str(exc)))

# Actual fd/scandir/stat implementation, only tiny private files and injected I/O.
original_os=m.os
for mode in ['healthy','removal','permission','EIO','cap','symlink']:
    root=FIX/('meter_'+mode);root.mkdir();q=root/'one';q.write_bytes(b'123')
    if mode=='symlink':(root/'link').symlink_to(q)
    class Entry:
        def __init__(self,e):self.e=e;self.name=e.name
        def stat(self,follow_symlinks):
            check(follow_symlinks is False)
            if self.name=='one':
                if mode=='removal':q.unlink()
                if mode=='permission':raise PermissionError(errno.EACCES,'INVENTED_PERMISSION')
                if mode=='EIO':raise OSError(errno.EIO,'INVENTED_IO')
            return self.e.stat(follow_symlinks=False)
    class Scan:
        def __init__(self,fd):self.real=os.scandir(fd)
        def __enter__(self):return iter(Entry(e) for e in self.real)
        def __exit__(self,*a):self.real.close()
    class OSProxy:
        def __getattr__(self,k):return getattr(original_os,k)
        def scandir(self,fd):return Scan(fd)
    m.os=OSProxy();m.SSD=root
    buf=io.StringIO()
    with contextlib.redirect_stdout(buf):
        value=case('meter_'+mode,lambda:m.global_namespace_gate({'SSD_reservation_bytes':2 if mode=='cap' else 10}),
                   PermissionError if mode=='permission' else OSError if mode=='EIO' else RuntimeError if mode=='cap' else None)
    if mode=='removal':check(value==0 and 'GLOBAL_METER_CONCURRENT_ENOENT_REMOVAL' in buf.getvalue())
    if mode in ['healthy','symlink']:check(value==3)
    if mode in ['permission','EIO']:check('GLOBAL_METER_CONCURRENT_ENOENT_REMOVAL' not in buf.getvalue())
m.os=original_os

# Execute the unchanged acquisition lifecycle with inert Popen and 15-byte body.
full=b'abcdeFGHIJKLMNO';prefix=full[:5]
md5=lambda b:hashlib.md5(b).hexdigest()
sha=lambda b:hashlib.sha256(b).hexdigest()
for mode in ['healthy206','wrong200','wrong_range','wrong_etag','wrong_full_sha','wrong_prefix']:
    root=FIX/('range_'+mode);root.mkdir();m.SSD=root;m.FOLDER=root/'new';m.PACKAGE=root/'package'
    for sub in ['raw','logs','receipts']:(m.FOLDER/sub).mkdir(parents=True,exist_ok=True)
    old=root/'preserved.partial';old.write_bytes(prefix if mode!='wrong_prefix' else b'zzzzz')
    member=dict(index=12,source_id='invented_source12',filename='invented.gz',body_path=str(m.FOLDER/'raw/invented.gz'),
        expected_bytes=len(full),expected_md5=md5(full),expected_sha256=sha(full),
        url='https://invented.invalid/body?generation=123',original_http_etag='"'+md5(full)+'"',
        expected_transport_headers={'etag':md5(full),'x-goog-generation':'123','x-goog-stored-content-length':str(len(full))})
    plan=dict(resume_source12=dict(original_partial_path=str(old),prefix_bytes=5,prefix_sha256=sha(prefix)),
        SSD_reservation_bytes=1<<20,internal_floor_bytes=3<<30,ssd_floor_bytes=5<<30,
        per_body_seconds_limit=7200,family_seconds_limit=96*3600,runtime_poll_seconds=2)
    calls=[]
    class Proc:
        pid=123456;returncode=0
        def poll(self):return 0
        def wait(self,timeout):return 0
    def popen(command,**kwargs):
        calls.append(command)
        check(command[:2]==['/usr/bin/curl','-q'])
        check(command[command.index('--continue-at')+1]=='5')
        check(command[command.index('--header')+1]=='If-Match: '+member['original_http_etag'])
        check('--retry' not in command)
        partial=Path(command[command.index('--output')+1]);check(partial.read_bytes()==prefix)
        with partial.open('ab') as f:f.write(full[5:] if mode!='wrong_full_sha' else b'XXXXXXXXXX')
        header=Path(command[command.index('--dump-header')+1])
        fields={'content-length':'10','content-range':'bytes 5-14/15',**member['expected_transport_headers']}
        if mode=='wrong_range':fields['content-range']='bytes 6-14/15'
        if mode=='wrong_etag':fields['etag']='wrong'
        header.write_text('HTTP/2 '+('200' if mode=='wrong200' else '206')+'\n'+''.join(k+': '+v+'\n' for k,v in fields.items()))
        return Proc()
    mon=SimpleNamespace(snapshot=lambda:dict(internal_free_bytes=100<<30,ssd_free_bytes=100<<30),
        final_limits=lambda *a:None,group_members=lambda pid:[],
        terminate_owned=lambda proc:dict(teardown_verified=True,remaining_group_members=[],signals=[]))
    m.subprocess=SimpleNamespace(Popen=popen,STDOUT=-2);m.physical_mount=lambda:None;m.assert_bindings=lambda *a:None
    class PopenOSProxy:
        def __getattr__(self,k):return getattr(original_os,k)
        def getpgid(self,pid):return pid
    m.os=PopenOSProxy();m.TERMINATION_REQUEST.clear()
    with contextlib.redirect_stdout(io.StringIO()):
        item=case('range_'+mode,lambda:m.acquire(plan,member,mon,'0'*64,time.monotonic()),None if mode=='healthy206' else RuntimeError)
    check(old.read_bytes()==(prefix if mode!='wrong_prefix' else b'zzzzz'))
    check(len(calls)==(0 if mode=='wrong_prefix' else 1))
    receipt=json.loads((m.FOLDER/'receipts/invented_source12.json').read_text())
    if mode=='healthy206':
        check(Path(member['body_path']).read_bytes()==full)
        check(receipt['actual_sha256']==sha(full) and receipt['resume_offset']==5)
    else:
        check(receipt['status']=='FAILED_PRESERVED_NO_AUTOMATIC_RETRY')
        check(not Path(member['body_path']).exists())
    m.os=original_os

# Focused mixed-origin gate: actual method, tiny source, no actual11 rehash.
root=FIX/'origin';root.mkdir();body=root/'body';body.write_bytes(full)
m.PACKAGE=P;m.CURL=Path('/usr/bin/curl')
header=root/'headers';header.write_text('HTTP/2 200\ncontent-length: 15\netag: invented\n')
member=dict(index=1,source_id='invented_reused',body_path=str(body),url='https://invented.invalid/original',expected_bytes=15,
    expected_md5=md5(full),expected_sha256=sha(full),expected_transport_headers={'etag':'invented'})
plan=dict(per_body_seconds_limit=7200,resume_source12={'prefix_bytes':5},receipt_origins={})
receipt=dict(member=member,status='EXACT_ORIGINAL_SOURCE_BODY_ACQUIRED',plan_sha256=m.ORIGINAL_PLAN_SHA,returncode=0,stop_reason=None,
    teardown=dict(teardown_verified=True,remaining_group_members=[]),resume_offset=0,post_cleanup_hash_resource_identity_gates_pass=True,
    actual_size=15,actual_md5=md5(full),actual_sha256=sha(full),final_seal_md5=md5(full),final_seal_sha256=sha(full),
    headers_sha256=sha(header.read_bytes()),observed_headers=m.header_fields(header),
    command=m.transfer_command(plan,member,header,str(body)+'.partial'))
for mode in ['healthy','wrong_origin_plan','cleanup_error','wrong_command','failed_status','failure_marker']:
    d=root/mode;d.mkdir();primary=d/'receipt';mirror=d/'mirror';r=copy.deepcopy(receipt)
    if mode=='wrong_origin_plan':r['plan_sha256']='0'*64
    if mode=='cleanup_error':r['teardown']['cleanup_error']='INVENTED'
    if mode=='wrong_command':r['command'].append('--retry')
    if mode=='failed_status':r['status']='FAILED_PRESERVED_NO_AUTOMATIC_RETRY'
    payload=json.dumps(r,indent=2)+'\n';primary.write_text(payload);mirror.write_text(payload);digest=sha(payload.encode())
    if mode=='failure_marker':Path(str(primary)+'.failure.json').write_text('{}')
    plan['receipt_origins'][member['source_id']]=dict(reused_v3=True,primary_path=str(primary),mirror_path=str(mirror),header_path=str(header),sha256=digest)
    case('reuse_'+mode,lambda:m.source_receipt_gate(plan,member,dict(path=str(primary),sha256=digest),'1'*64),None if mode=='healthy' else RuntimeError)

# Unchanged ownership/Terminal2 helpers are inherited, not rerun.
a=ast.parse((P/'scripts/90_acquire_validation_raw_sources_v3.py').read_text());b=ast.parse(CODE.read_text())
fa={n.name:ast.dump(n,include_attributes=False) for n in a.body if isinstance(n,ast.FunctionDef)}
fb={n.name:ast.dump(n,include_attributes=False) for n in b.body if isinstance(n,ast.FunctionDef)}
names=['utc','safe_print','catchable_termination','hashes','write_new','family_lock','load_monitor','physical_mount','assert_bindings','header_fields','safe_state','cleanup_proof','retain_lock_until_gone','regular']
check(all(fa[n]==fb[n] for n in names))
check('TerminalCommit' in CODE.read_text())
files={str(q):sha(q.read_bytes()) for q in FIX.rglob('*') if q.is_file() and not q.is_symlink()}
check(sum(q.stat().st_size for q in FIX.rglob('*') if q.is_file() and not q.is_symlink())<128<<10)
check(time.monotonic()-started<30 and resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<128<<20)
r=dict(schema='validation13_raw_v4_author_focused_delta_controls',control_count=len(results),assertion_count=checks,
    controls=results,unchanged_helper_AST=names,code_sha256=sha(CODE.read_bytes()),fixture_file_sha256=files,
    no_production_body_reads=True,no_network_or_actual_workers_or_locks=True,no_old_suites_rerun=True,
    elapsed_seconds=time.monotonic()-started,max_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
with OUT.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
print(json.dumps(dict(controls=len(results),assertions=checks,receipt_sha256=sha(OUT.read_bytes()),elapsed=r['elapsed_seconds'])))

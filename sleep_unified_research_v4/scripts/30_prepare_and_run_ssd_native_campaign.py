#!/usr/bin/env python3
"""Prepare current hashes and monitor unchanged native stages in a new SSD namespace."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / 'sleep_unified_research_v4'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
SUPPORT = SSD / 'native' / 'original_runner_support_v4_2'
PLAN_PATH = PACKAGE / 'manifests/ssd_native_execution_plan_v4_3.json'
INTERNAL_FLOOR = 3 * 1024**3
SSD_FLOOR = 5 * 1024**3
RSS_LIMIT = 2 * 1024**3
STAGE_SECONDS = 36 * 3600
OUTPUT_LIMIT = 4 * 1024**3

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()

def rows(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def utc():
    return datetime.now(timezone.utc).isoformat()

def snapshot():
    return {'recorded_utc': utc(),
            'internal_free_bytes': shutil.disk_usage('/System/Volumes/Data').free,
            'ssd_free_bytes': shutil.disk_usage(SSD).free,
            'internal_minimum_bytes': INTERNAL_FLOOR, 'ssd_minimum_bytes': SSD_FLOOR}

def module():
    path = ROOT / 'sleep_unified_research_v1/scripts/04_native_reproduction_runner.py'
    assert sha(path) == 'fa5ca8528caa1037cca55d5234e7e7721582abce479c95949bb81255e3428fa9'
    spec = importlib.util.spec_from_file_location('sleep_native_v1_prepare', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def prepare():
    if PLAN_PATH.exists():
        raise RuntimeError('FROZEN_PREPARATION_ALREADY_EXISTS')
    mod = module()
    jobs = mod.jobs()
    original_plan_path = ROOT / 'sleep_unified_research_v1/manifests/native_reproduction_jobs_v1.json'
    original = json.loads(original_plan_path.read_text())
    assert jobs == original['jobs'] and len(jobs) == 190
    inputs = sorted({source for job in jobs for source in job['inputs']})
    admitted = {r['path']:r['actual_sha256'] for r in rows(ROOT / 'sleep_unified_research_v1/tables/native_input_hash_checks.tsv') if r['kind']=='core_munged'}
    admitted.update({r['path']:r['actual_sha256'] for r in rows(ROOT / 'sleep_unified_research_v1/tables/archived_extension_recovery.tsv') if r['status']=='EXACT_RECEIPT_HASH_RECOVERED'})
    assert len(inputs)==158 and all(path in admitted for path in inputs)
    SUPPORT.mkdir(parents=True, exist_ok=False)
    support_hashes = {}
    for name in ['FROZEN_NEW_ANALYSIS_PROTOCOL.md','scripts/native_ldsc_capture.py','tables/native_input_hash_checks.tsv','tables/archived_extension_recovery.tsv']:
        source = ROOT / 'sleep_unified_research_v1' / name
        target = SUPPORT / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        assert sha(source)==sha(target)
        support_hashes[str(target)] = sha(target)
    for directory in ('manifests','logs','native'):
        (SUPPORT / directory).mkdir(exist_ok=True)
    checks = []
    for number, source in enumerate(inputs,1):
        if snapshot()['internal_free_bytes'] < 128 * 1024**2:
            raise RuntimeError('STDLIB_HASH_EMERGENCY_FLOOR')
        actual = sha(source)
        checks.append({'path':source,'bytes':Path(source).stat().st_size,
                       'expected_sha256':admitted[source],'actual_sha256':actual,
                       'match':actual==admitted[source]})
        assert checks[-1]['match'], source
        if number % 25 == 0:
            print(json.dumps({'processed_inputs_verified':number,'total':len(inputs)}),flush=True)
    deps = mod.execution_dependencies()
    env_command = [str(mod.PYTHON),'-B','-c','import json,sys,numpy,pandas,scipy;print(json.dumps(dict(python=sys.version,numpy=numpy.__version__,pandas=pandas.__version__,scipy=scipy.__version__)))']
    versions = subprocess.run(env_command,capture_output=True,text=True,check=True,
                              env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','TMPDIR':str(SSD/'tmp'),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
    environment = json.loads(versions.stdout)
    assert environment['python'].startswith('3.9.23')
    assert (environment['numpy'],environment['pandas'],environment['scipy'])==('1.21.5','1.3.3','1.7.3')
    plan = {'prepared_utc':utc(),'jobs':jobs,'original_plan_sha256':sha(original_plan_path),
            'original_runner_sha256':sha(ROOT/'sleep_unified_research_v1/scripts/04_native_reproduction_runner.py'),
            'launcher_sha256':sha(Path(__file__)),
            'invoker_sha256':sha(PACKAGE/'scripts/31_invoke_original_native_runner.py'),
            'ssd_support_package':str(SUPPORT),'support_file_sha256':support_hashes,
            'dependencies_sha256':deps,'inputs_verified':checks,'environment':environment,
            'worker_count':1,'blas_threads':1,'internal_floor_bytes':INTERNAL_FLOOR,
            'ssd_floor_bytes':SSD_FLOOR,'maximum_observed_aggregate_worker_rss_bytes':RSS_LIMIT,
            'stage_seconds_limit':STAGE_SECONDS,'new_native_output_byte_limit':OUTPUT_LIMIT,
            'runtime_poll_seconds':2,'operational_change':'Explicit PACKAGE output/support relocation only; estimator ROOT, job identities, options, inputs, original runner and resource floor unchanged.',
            'no_new_primary_or_independent_test_admitted':True,'resource_preflight':snapshot()}
    PLAN_PATH.write_text(json.dumps(plan,indent=2)+'\n')
    (SSD/'manifests'/PLAN_PATH.name).write_bytes(PLAN_PATH.read_bytes())
    print(json.dumps({'jobs_prepared':len(jobs),'current_inputs_verified':len(inputs),
                      'dependencies':len(deps),'resource_preflight':plan['resource_preflight'],
                      'native_jobs_launched':0}),flush=True)

def owned_rss(pid):
    result=subprocess.run(['ps','-axo','pid=,ppid=,rss='],capture_output=True,text=True,check=True)
    processes=[tuple(map(int,line.split())) for line in result.stdout.splitlines() if len(line.split())==3]
    owned={pid}
    while True:
        more={p for p,parent,rss in processes if parent in owned}
        if more <= owned:break
        owned |= more
    return sum(rss*1024 for p,parent,rss in processes if p in owned),sorted(owned)

def group_members(group_id):
    result=subprocess.run(['ps','-axo','pid=,pgid=,stat='],capture_output=True,text=True,check=True)
    return [{'pid':int(parts[0]),'state':parts[2]} for line in result.stdout.splitlines()
            if len(parts:=line.split())==3 and int(parts[1])==group_id]

def terminate_owned(proc):
    # Group persists after its immediate parent exits; track it until empty.
    initial=group_members(proc.pid);signals=[]
    if initial:
        try:os.killpg(proc.pid,signal.SIGTERM);signals.append('SIGTERM')
        except ProcessLookupError:pass
        deadline=time.monotonic()+10
        while time.monotonic()<deadline:
            proc.poll()
            if not group_members(proc.pid):break
            time.sleep(.25)
        if group_members(proc.pid):
            try:os.killpg(proc.pid,signal.SIGKILL);signals.append('SIGKILL')
            except ProcessLookupError:pass
            deadline=time.monotonic()+10
            while time.monotonic()<deadline:
                proc.poll()
                if not group_members(proc.pid):break
                time.sleep(.25)
    proc.wait(timeout=10)
    remaining=group_members(proc.pid)
    if remaining:raise RuntimeError('OWNED_PROCESS_GROUP_REMAINS_AFTER_TEARDOWN')
    return {'initial_group_members':initial,'signals':signals,'remaining_group_members':remaining}

def final_limits(state,rss,output,elapsed):
    if state['internal_free_bytes'] < INTERNAL_FLOOR:return 'INTERNAL_FULL_NATIVE_FLOOR_REACHED'
    if state['ssd_free_bytes'] < SSD_FLOOR:return 'SSD_FULL_NATIVE_FLOOR_REACHED'
    if rss > RSS_LIMIT:return 'WORKER_RSS_LIMIT'
    if output > OUTPUT_LIMIT:return 'NATIVE_OUTPUT_LIMIT'
    if elapsed > STAGE_SECONDS:return 'STAGE_DEADLINE'
    return None

def execute(stage):
    plan=json.loads(PLAN_PATH.read_text())
    assert plan['launcher_sha256']==sha(Path(__file__))
    assert plan['invoker_sha256']==sha(PACKAGE/'scripts/31_invoke_original_native_runner.py')
    state=snapshot()
    receipt_path=PACKAGE/'logs'/f'{stage}_native_monitor_receipt_v4.json'
    if receipt_path.exists():raise RuntimeError('PRIOR_STAGE_MONITOR_REQUIRES_REVIEW')
    if state['internal_free_bytes'] < INTERNAL_FLOOR or state['ssd_free_bytes'] < SSD_FLOOR:
        (PACKAGE/'logs'/f'{stage}_launch_guard_v4.json').write_text(json.dumps({**state,'native_jobs_launched':0,'resource_gate_pass':False},indent=2)+'\n')
        raise SystemExit('FULL_NATIVE_RESOURCE_GATE_FAILED_NO_WORKER_LAUNCHED')
    assert module().execution_dependencies()==plan['dependencies_sha256']
    command=[sys.executable,'-B',str(PACKAGE/'scripts/31_invoke_original_native_runner.py'),stage]
    env={**os.environ,'TMPDIR':str(SSD/'tmp'),'XDG_CACHE_HOME':str(SSD/'cache'),
         'MPLCONFIGDIR':str(SSD/'cache/matplotlib'),'PYTHONDONTWRITEBYTECODE':'1',
         'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
    stdout=SSD/'logs'/f'{stage}_original_runner_stdout_v4.log'
    started=time.monotonic();proc=None;stop_reason=None;samples=[];peak=0;teardown=None
    receipt={'stage':stage,'command':command,'started_utc':utc(),'plan_sha256':sha(PLAN_PATH),
             'resource_preflight':state,'env_paths':{k:env[k] for k in ['TMPDIR','XDG_CACHE_HOME','MPLCONFIGDIR']}}
    try:
        with stdout.open('x') as log:
            proc=subprocess.Popen(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            assert os.getpgid(proc.pid)==proc.pid,'NEW_OWNED_SESSION_IDENTITY_FAILED'
            receipt['worker_pid']=proc.pid
            while proc.poll() is None:
                current=snapshot();rss,pids=owned_rss(proc.pid);peak=max(peak,rss)
                output=sum(p.stat().st_size for p in (SUPPORT/'native').rglob('*') if p.is_file())
                sample={**current,'aggregate_owned_rss_bytes':rss,'owned_pids':pids,'native_output_bytes':output}
                samples.append(sample)
                stop_reason=final_limits(current,rss,output,time.monotonic()-started)
                if stop_reason:
                    teardown=terminate_owned(proc);break
                if len(samples)%30==0:
                    print(json.dumps({'stage':stage,'pid':proc.pid,'elapsed_seconds':round(time.monotonic()-started,1),'rss_bytes':rss,'internal_free_bytes':current['internal_free_bytes']}),flush=True)
                time.sleep(2)
    except BaseException as exc:
        stop_reason=stop_reason or type(exc).__name__+': '+str(exc)
        raise
    finally:
        if proc is not None:
            try:
                teardown=teardown or terminate_owned(proc)
                if group_members(proc.pid):raise RuntimeError('OWNED_GROUP_REAPPEARED')
                if teardown['initial_group_members'] and stop_reason is None:
                    stop_reason='UNEXPECTED_DESCENDANTS_AFTER_INVOKER_EXIT'
            except BaseException as cleanup_error:
                stop_reason='PROCESS_TEARDOWN_FAILURE: '+str(cleanup_error)
                teardown={'cleanup_error':str(cleanup_error),'teardown_verified':False}
        final_state=snapshot()
        final_output=sum(p.stat().st_size for p in (SUPPORT/'native').rglob('*') if p.is_file())
        final_elapsed=time.monotonic()-started
        stop_reason=stop_reason or final_limits(final_state,peak,final_output,final_elapsed)
        receipt.update({'completed_utc':utc(),'elapsed_seconds':time.monotonic()-started,
                        'returncode':proc.returncode if proc else None,'stop_reason':stop_reason,
                        'process_group_teardown':teardown,'final_resource_snapshot':final_state,
                        'final_native_output_bytes':final_output,
                        'peak_observed_owned_rss_bytes':peak,'samples':samples,
                        'stdout_path':str(stdout),'stdout_sha256':sha(stdout) if stdout.exists() else None,
                        'native_receipt_paths':[str(p) for p in (SUPPORT/'native').rglob('*.execution_receipt.json')],
                        'scientific_completion_requires_per_job_and_numerical_review':True})
        receipt_path.write_text(json.dumps(receipt,indent=2)+'\n')
        (SSD/'logs'/receipt_path.name).write_bytes(receipt_path.read_bytes())
    if receipt['returncode']!=0 or stop_reason:raise SystemExit('NATIVE_STAGE_STOP_PRESERVED_REQUIRES_REVIEW')
    print(json.dumps({'stage':stage,'native_receipts':len(receipt['native_receipt_paths']),'returncode':receipt['returncode']}),flush=True)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--prepare',action='store_true')
    parser.add_argument('--stage',choices=['core','extension','validation']);parser.add_argument('--execute',action='store_true')
    args=parser.parse_args()
    if args.prepare:prepare()
    elif args.execute and args.stage:execute(args.stage)
    else:parser.error('Use --prepare or --stage STAGE --execute')

if __name__=='__main__':main()

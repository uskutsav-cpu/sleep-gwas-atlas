#!/usr/bin/env python3
"""Independent bounded descriptor/admission software controls; no GWAS reads."""
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[2]
P = ROOT/'sleep_unified_research_v4'
sys.path.insert(0, str(P/'scripts'))
SOURCE = P/'scripts/58_run_finngen_feasibility_stage.py'
RECEIPT = P/'reviews/independent_finngen_launcher_proxy_controls_receipt_v4.json'
TEMP_PARENT = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/reviewer_software_controls')

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    started = time.monotonic()
    before = sha(SOURCE)
    spec = importlib.util.spec_from_file_location('_independent_finngen_adapter', SOURCE)
    candidate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(candidate)
    original_popen = subprocess.Popen
    TEMP_PARENT.mkdir(exist_ok=True)
    checks = []
    proc = None
    with tempfile.TemporaryDirectory(prefix='independent-finn-proxy-',dir=TEMP_PARENT) as td:
        directory = Path(td)
        lock = directory/'temporary.lock'
        fd = os.open(lock, os.O_RDWR | os.O_CREAT, 0o600)
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        expected = os.fstat(fd)
        proxy = candidate.InheritedMutexSubprocess(fd)
        child = ('import fcntl,json,os,sys; f=int(sys.argv[1]); s=os.fstat(f); '
                 'fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB); '
                 'print(json.dumps([s.st_dev,s.st_ino]),flush=True); sys.stdin.readline()')
        contender = ('import fcntl,sys; f=open(sys.argv[1],"a+"); '
                     '\ntry: fcntl.flock(f.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB); print("FREE")'
                     '\nexcept BlockingIOError: print("BUSY")')
        try:
            proc = proxy.Popen([sys.executable,'-B','-c',child,str(fd)],stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,start_new_session=True)
            # Read readiness with a bounded poll before touching the pipe.
            import select
            if not select.select([proc.stdout], [], [], 3)[0]:
                raise RuntimeError('HELPER_READINESS_TIMEOUT')
            identity = json.loads(proc.stdout.readline())
            assert identity == [expected.st_dev, expected.st_ino]
            checks.append('inherited descriptor exact device/inode and flock reacquisition')
            os.close(fd); fd = None
            busy = subprocess.run([sys.executable,'-B','-c',contender,str(lock)],capture_output=True,text=True,timeout=3,check=True)
            assert busy.stdout.strip() == 'BUSY'
            checks.append('child alone retains mutex after parent descriptor closes')
            proc.communicate('exit\n',timeout=3)
            assert proc.returncode == 0
            free = subprocess.run([sys.executable,'-B','-c',contender,str(lock)],capture_output=True,text=True,timeout=3,check=True)
            assert free.stdout.strip() == 'FREE'
            checks.append('mutex available only after inherited helper exits')
            try:
                proxy.Popen([sys.executable,'-c','pass'],pass_fds=())
            except RuntimeError as error:
                assert str(error) == 'UNEXPECTED_WORKER_DESCRIPTOR_OVERRIDE'
            else:
                raise AssertionError('descriptor override accepted')
            checks.append('explicit pass_fds override rejected before child creation')
            assert subprocess.Popen is original_popen
            checks.append('global subprocess module unmodified')
        finally:
            if fd is not None:os.close(fd)
            if proc is not None and proc.poll() is None:
                os.killpg(proc.pid, signal.SIGTERM)
                try:proc.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL);proc.wait(timeout=2)
        # Admission-only fixture uses synthetic metadata, never worker execution.
        original = json.loads((P/'manifests/finngen_preprocessing_plan_v4.json').read_text())
        base = directory/'synthetic-baseline.json';base.write_text('{}\n')
        admission = directory/'synthetic-admission.json'
        plan = dict(original, dependencies_sha256={}, baseline_gate_plan=str(base),
                    baseline_gate_plan_sha256=sha(base), admission=str(admission))
        plan_path = directory/'synthetic-plan.json';plan_path.write_text(json.dumps(plan))
        plan_hash = sha(plan_path)
        admission.write_text(json.dumps(dict(execution_admitted=True,plan_sha256=plan_hash,
            scope=plan['scope'],independent_binding_review_pass=True,resource_plan_review_pass=True,
            independent_review_sha256={})))
        candidate.require_admission(plan_path, plan_hash)
        empty_map_accepted = True
    assert sha(SOURCE) == before
    receipt = dict(status='PASS_PROXY_CONTROLS_WITH_CONFIRMED_EMPTY_REVIEW_MAP_ADMISSION_FAULT',
        executor_sha256=before,checker_sha256=sha(__file__),checks=checks,
        empty_independent_review_hash_map_accepted=empty_map_accepted,
        actual_helpers=3,actual_GWAS_or_preprocessing_or_native_workers=0,
        real_study_lock_used=False,GWAS_body_reads=0,source_code_unchanged=True,
        elapsed_seconds=time.monotonic()-started,
        scope='Temporary independent lock and tiny stdlib helpers only; synthetic admission metadata.')
    with RECEIPT.open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
    print(json.dumps(receipt))

if __name__ == '__main__':main()

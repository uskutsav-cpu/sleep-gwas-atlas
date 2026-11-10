"""Validate actual process state and prevent duplicate finite-comparison jobs."""
import contextlib,json,os,subprocess,time
from pathlib import Path


def command(pid):
    result=subprocess.run(['ps','-p',str(pid),'-o','command='],capture_output=True,text=True)
    return result.stdout.strip() if result.returncode==0 else ''


def active_fit(script,recipe,source):
    rows=subprocess.run(['ps','-axo','pid,command'],capture_output=True,text=True,check=True).stdout.splitlines()
    result=[]
    for row in rows:
        fields=row.strip().split(None,1)
        if len(fields)!=2:continue
        pid,cmd=fields
        if script in cmd and f' --recipe {recipe} ' in ' '+cmd+' ' and f' --held-study {source} ' in ' '+cmd+' ':result.append(int(pid))
    assert len(result)<=1,('Duplicate training processes',result)
    return result[0] if result else None


@contextlib.contextmanager
def driver_lock(path):
    path=Path(path);needle='run_resumable_comparison.py'
    if path.exists():
        previous=json.loads(path.read_text());cmd=command(previous['PID'])
        if needle in cmd:raise SystemExit(f'Existing driver confirmed live: {previous["PID"]}; follow that process.')
        path.unlink()
    descriptor=os.open(path,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
    with os.fdopen(descriptor,'w') as handle:json.dump({'PID':os.getpid(),'script':needle},handle)
    try:yield
    finally:
        if path.exists() and json.loads(path.read_text())['PID']==os.getpid():path.unlink()

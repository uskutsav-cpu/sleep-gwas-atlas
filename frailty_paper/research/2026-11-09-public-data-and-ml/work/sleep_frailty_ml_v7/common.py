"""Separate V7 workspace; previous scientific locks and models stay immutable."""
from pathlib import Path
import os,json,csv,hashlib,time,ssl,urllib.request
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
WORK=Path(__file__).resolve().parent
DATA=Path('/Volumes/Extreme SSD/codex-sleep-frailty-work/sleep_frailty_ml_v7')
OUT=DATA/'outputs';SOURCES=DATA/'sources';MODELS=DATA/'models'
V6=DATA.parent/'sleep_frailty_ml_v6';V6OUT=V6/'outputs';OLD=ROOT/'work/sleep_frailty_ml'
SSL=ssl.create_default_context(cafile='/etc/ssl/cert.pem')
for p in [DATA,OUT,SOURCES,MODELS]:p.mkdir(parents=True,exist_ok=True)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def write_json(path,value):
    p=Path(path);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(p)
def fetch(url,path,limit=100_000_000,magic=None):
    path=Path(path);temp=path.with_suffix(path.suffix+'.download')
    with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'SleepFrailtyResearch/1.0'}),context=SSL,timeout=60) as r,temp.open('wb') as f:
        n=0
        while True:
            b=r.read(1024*1024)
            if not b:break
            n+=len(b);assert n<=limit,'Source exceeds declared download bound';f.write(b)
    if magic is not None:
        with temp.open('rb') as f:assert f.read(len(magic))==magic,'Unexpected source format; not a valid data artifact'
    temp.replace(path);return {'url':url,'path':str(path),'bytes':path.stat().st_size,'sha256':sha(path)}
def environment():
    env=os.environ.copy();alias=Path('/tmp/codex_sf_v7_runtime');folder=DATA/'runtime_tmp';folder.mkdir(exist_ok=True)
    if alias.exists() or alias.is_symlink():assert alias.is_symlink() and alias.resolve()==folder.resolve()
    else:alias.symlink_to(folder,target_is_directory=True)
    env['TMPDIR']=str(alias);env['PYTHONDONTWRITEBYTECODE']='1'
    for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:env[k]='2'
    return env
def assert_previous_intact():
    lock=json.loads((V6OUT/'final_model_lock.json').read_text())
    for p,h in {**lock['files'],**lock['model_checkpoints']}.items():assert sha(p)==h,p
    assert sha(ROOT/'outputs/sleep_frailty_ai_v3/final_model_registry.json')==lock['V3_registry_sha256']
    return sha(V6OUT/'final_model_lock.json')

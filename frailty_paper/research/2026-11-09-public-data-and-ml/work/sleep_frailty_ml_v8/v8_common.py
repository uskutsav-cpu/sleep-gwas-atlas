"""New measured-effect experiment. V7 and prior evidence remain immutable."""
import csv, gzip, hashlib, json, os, ssl, sys, urllib.request
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
BASE = Path('/Volumes/Extreme SSD/codex-sleep-frailty-work')
DATA = BASE/'sleep_frailty_ml_v8'
OUT = DATA/'outputs'
SOURCES = DATA/'sources'
MODELS = DATA/'models'
V7 = BASE/'sleep_frailty_ml_v7'
V6 = BASE/'sleep_frailty_ml_v6'
for p in [DATA, OUT, SOURCES, MODELS]:
    p.mkdir(parents=True, exist_ok=True)
alias = ROOT/'outputs/sleep_frailty_ai_v8'
if not alias.exists():
    alias.symlink_to(OUT, target_is_directory=True)
assert alias.resolve() == OUT.resolve()
SSL = ssl.create_default_context(cafile='/etc/ssl/cert.pem')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value, indent=2)+'\n')
    tmp.replace(path)


def fetch(url, path, limit, magic=None):
    path = Path(path)
    if not path.exists():
        tmp = path.with_suffix(path.suffix+'.download')
        with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'SleepFrailtyResearch/1.0'}), context=SSL, timeout=45) as r, tmp.open('wb') as f:
            n = 0
            while True:
                b = r.read(1024*1024)
                if not b:
                    break
                n += len(b)
                assert n <= limit, 'Download exceeds declared bound'
                f.write(b)
        if magic:
            with tmp.open('rb') as f:
                assert f.read(len(magic)) == magic
        tmp.replace(path)
    assert path.stat().st_size <= limit
    return {'url': url, 'path': str(path), 'bytes': path.stat().st_size, 'sha256': sha(path)}


def rc(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def csv_rows(path):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt') as f:
        yield from csv.DictReader(f)


def write_csv(path, rows):
    with Path(path).open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def assert_previous_intact():
    lock = json.loads((V6/'outputs/final_model_lock.json').read_text())
    for p, h in {**lock['files'], **lock['model_checkpoints']}.items():
        assert sha(p) == h, p
    assert sha(ROOT/'outputs/sleep_frailty_ai_v3/final_model_registry.json') == lock['V3_registry_sha256']
    receipt = json.loads((V7/'outputs/scientific_report_refinement_receipt.json').read_text())
    for name, h in receipt['frozen_files_unchanged'].items():
        assert sha(V7/'outputs'/name) == h, name


def environment():
    env = os.environ.copy()
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS']:
        env[k] = '2'
    return env

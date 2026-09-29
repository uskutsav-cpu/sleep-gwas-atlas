#!/usr/bin/env python3
"""Freeze and run isolated postmortem of v2 insomnia worker 2 only."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'brain6/results/confirmatory_source_rescue_20260927/other_traits'
PROTOCOL = BASE / 'INSOMNIA_NATIVE_N_WORKER2_POSTMORTEM_V3.md'
CONFIG = ROOT / 'brain6/config/confirmatory_source_rescue_20260927/insomnia_native_n_worker2_postmortem_v3.json'
RUNNER = ROOT / 'brain6/scripts/diagnose_insomnia_native_n_worker2_v3.R'
RSCRIPT = Path('/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin/Rscript')
V2 = Path('/Volumes/Extreme SSD/brain6-work/confirmatory-source-rescue-20260927/insomnia-native-n-pilot-v2')
V3 = Path('/Volumes/Extreme SSD/brain6-work/confirmatory-source-rescue-20260927/insomnia-native-n-worker2-postmortem-v3')
V2_CONFIG = ROOT / 'brain6/config/confirmatory_source_rescue_20260927/insomnia_native_n_pilot_v2.json'
EXPECTED = {
    V2_CONFIG: '893b0b8251805a421355e3141b658c1997bed4783daf7990a378abeb32796401',
    V2 / 'materialization.receipt.json': '24ce34356859f1a3f58390c1dc9193c6849cabe2944b1ac947b08d23f53e4f83',
    V2 / 'results/launch.json': 'b5a3b7e0d60f43a38611150956e0fa6b302a40197885a88c16a7ea3d30039060',
    V2 / 'results/exit.json': '4e3c578c526353d3a0ef8108b9e2915f8e2a47c599e0b30f0358b7a7bcfe1aa1',
    V2 / 'results/worker_2.log': '0bd8e5f5f28ca35e58d6fdb5764b0bee54e0d0f3d93d0985db3fa57ba8336830',
}

def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()

def check(path: Path, expected: str) -> None:
    if sha(path) != expected:
        raise ValueError(f'Frozen hash mismatch: {path}')

def save(path: Path, obj: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(obj, stream, indent=2, sort_keys=True)
        stream.write('\n')

def validate_v2() -> tuple[dict, dict]:
    for path, digest in EXPECTED.items():
        check(path, digest)
    if (V2 / 'results/worker_2.tsv').exists() or (V2 / 'results/worker_2.summary.json').exists():
        raise ValueError('V2 worker 2 output unexpectedly exists')
    exit_data = json.loads((V2 / 'results/exit.json').read_text())
    if exit_data['exit_codes'] != [0, 1, 0, 0]:
        raise ValueError('V2 failure status changed')
    worker = json.loads((V2 / 'worker_2.config.json').read_text())
    if worker['worker_id'] != 2 or len(worker['locus_ids']) != 22:
        raise ValueError('V2 worker 2 locus identity changed')
    receipt = json.loads((V2 / 'materialization.receipt.json').read_text())
    for loc in worker['locus_ids']:
        record = receipt['records'][loc]
        check(V2 / 'inputs' / f'locus_{loc}' / 'insomnia.sumstats.tsv.gz', record['sumstats_sha256'])
        check(V2 / 'inputs' / f'locus_{loc}' / 'input_info.tsv', record['input_info_sha256'])
    check(V2 / 'selected.loci', receipt['selected_loci_sha256'])
    return worker, receipt

def freeze() -> None:
    if CONFIG.exists() or V3.exists():
        raise FileExistsError('Postmortem config or output root already exists')
    worker, _ = validate_v2()
    v3 = dict(worker)
    v3.update({'analysis_id': 'brain6_insomnia_native_n_worker2_postmortem_v3',
               'scope': 'POSTMORTEM_ONLY_NO_ADVANCEMENT',
               'output_tsv': str(V3 / 'worker_2.diagnostic.tsv'),
               'summary_json': str(V3 / 'worker_2.diagnostic.summary.json'),
               'v2_worker_config_sha256': sha(V2 / 'worker_2.config.json'),
               'v2_materialization_receipt_sha256': EXPECTED[V2 / 'materialization.receipt.json'],
               'v2_exit_receipt_sha256': EXPECTED[V2 / 'results/exit.json'],
               'protocol_sha256': sha(PROTOCOL), 'r_runner_sha256': sha(RUNNER),
               'python_runner_sha256': sha(Path(__file__)), 'rscript_sha256': sha(RSCRIPT)})
    save(CONFIG, v3)
    print(json.dumps({'status': 'FROZEN_POSTMORTEM_ONLY', 'config_sha256': sha(CONFIG),
                      'worker': 2, 'loci': 22, 'output_root_unused': not V3.exists()}, sort_keys=True))

def run() -> None:
    cfg = json.loads(CONFIG.read_text())
    worker, _ = validate_v2()
    for key, path in [('protocol_sha256', PROTOCOL), ('r_runner_sha256', RUNNER),
                      ('python_runner_sha256', Path(__file__)), ('rscript_sha256', RSCRIPT)]:
        check(path, cfg[key])
    if cfg['v2_worker_config_sha256'] != sha(V2 / 'worker_2.config.json') or cfg['locus_ids'] != worker['locus_ids']:
        raise ValueError('Postmortem worker identity changed')
    if cfg['analysis_id'] != 'brain6_insomnia_native_n_worker2_postmortem_v3' or cfg['scope'] != 'POSTMORTEM_ONLY_NO_ADVANCEMENT':
        raise ValueError('Postmortem scope changed')
    if V3.exists():
        raise FileExistsError('Exclusive v3 output root exists')
    V3.mkdir(parents=True, exist_ok=False)
    launch = {'analysis_id': cfg['analysis_id'], 'scope': cfg['scope'],
              'launched_utc': datetime.now(timezone.utc).isoformat(),
              'config_sha256': sha(CONFIG), 'v2_exit_receipt_sha256': cfg['v2_exit_receipt_sha256'],
              'worker_id': 2, 'expected_loci': 22}
    save(V3 / 'launch.json', launch)
    env = dict(os.environ)
    env.update({'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1',
                'VECLIB_MAXIMUM_THREADS': '1', 'MKL_NUM_THREADS': '1'})
    log = V3 / 'worker_2.diagnostic.log'
    with log.open('xb') as stream:
        code = subprocess.run([str(RSCRIPT), str(RUNNER), str(CONFIG)], cwd=ROOT,
                              stdout=stream, stderr=subprocess.STDOUT, env=env, check=False).returncode
    results = V3 / 'worker_2.diagnostic.tsv'
    summary = V3 / 'worker_2.diagnostic.summary.json'
    if not results.exists() or not summary.exists():
        raise RuntimeError(f'Postmortem did not persist 22 diagnostic rows; R exit={code}')
    import csv
    with results.open(newline='') as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    if len(rows) != 22 or {r['locus_id'] for r in rows} != set(cfg['locus_ids']):
        raise ValueError('Postmortem output identity mismatch')
    failed = [{'locus_id': r['locus_id'], 'reason': r['reason']} for r in rows if r['status'] == 'FAILED']
    receipt = {'analysis_id': cfg['analysis_id'], 'scope': cfg['scope'],
               'status': 'REPRODUCED_FAILURE' if failed else 'FAILURE_NOT_REPRODUCED_V2_STILL_TERMINAL',
               'r_exit_code': code, 'rows': len(rows), 'failed': failed,
               'v2_advancement_rule_reopened': False, 'config_sha256': sha(CONFIG),
               'launch_sha256': sha(V3 / 'launch.json'), 'log_sha256': sha(log),
               'rows_sha256': sha(results), 'summary_sha256': sha(summary)}
    save(V3 / 'diagnostic.receipt.json', receipt)
    print(json.dumps({'status': receipt['status'], 'failed': failed,
                      'receipt_sha256': sha(V3 / 'diagnostic.receipt.json')}, sort_keys=True))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('freeze', 'run'))
    args = parser.parse_args()
    {'freeze': freeze, 'run': run}[args.action]()

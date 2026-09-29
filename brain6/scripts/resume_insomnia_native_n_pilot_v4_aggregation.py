#!/usr/bin/env python3
"""Complete v4's unchanged v2 decision rule without rerunning any worker."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / 'brain6/config/confirmatory_source_rescue_20260927/insomnia_native_n_pilot_v4.json'
V2_CODE = ROOT / 'brain6/scripts/run_insomnia_native_n_pilot_v2.py'
V2_ROOT = Path('/Volumes/Extreme SSD/brain6-work/confirmatory-source-rescue-20260927/insomnia-native-n-pilot-v2')
V4_ROOT = Path('/Volumes/Extreme SSD/brain6-work/confirmatory-source-rescue-20260927/insomnia-native-n-pilot-v4')
OUT = V4_ROOT / 'aggregation_repair_v1'


def main() -> None:
    cfg = json.loads(CONFIG.read_text())
    if OUT.exists():
        raise FileExistsError('Aggregation repair already exists; refuse overwrite')
    spec = importlib.util.spec_from_file_location('insomnia_v2_aggregator', V2_CODE)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    sha = mod.sha
    if sha(V2_CODE) != cfg['v2_python_runner_sha256'] or sha(CONFIG) != '9cd1f0a18a17e370e68e1e6170f8c9edda1a7b7ad86e5887505bd46e0d14d47c':
        raise ValueError('Frozen v4 config or original aggregator changed')
    if sha(V2_ROOT / 'materialization.receipt.json') != cfg['materialization_receipt_sha256']:
        raise ValueError('Original materialization receipt changed')
    rule = json.loads((ROOT / cfg['rule_relative']).read_text())
    if sha(ROOT / cfg['rule_relative']) != cfg['rule_sha256'] or rule['advance_if_any'] != cfg['advance_if_any']:
        raise ValueError('Frozen advancement rule changed')
    exit_info = json.loads((V4_ROOT / 'results' / 'exit.json').read_text())
    if exit_info['worker_2_exit_code'] != 0:
        raise ValueError('V4 worker 2 did not complete successfully')
    for suffix, key in (('.log', 'worker_2_log_sha256'), ('.summary.json', 'worker_2_summary_sha256'), ('.tsv', 'worker_2_tsv_sha256')):
        if sha(V4_ROOT / 'results' / f'worker_2{suffix}') != exit_info[key]:
            raise ValueError(f'V4 worker 2 {suffix} receipt mismatch')
    summary = json.loads((V4_ROOT / 'results' / 'worker_2.summary.json').read_text())
    if summary['rows'] != 22 or summary['failed'] != 0:
        raise ValueError('V4 worker 2 incomplete or failed')
    workers = []
    for worker in range(1, 5):
        path = V2_ROOT / f'worker_{worker}.config.json'
        if sha(path) != cfg['v2_worker_configs_sha256'][str(worker)]:
            raise ValueError(f'Original worker {worker} partition changed')
        source_root = V4_ROOT if worker == 2 else V2_ROOT
        if worker != 2:
            for suffix in ('.tsv', '.summary.json', '.log'):
                if sha(source_root / 'results' / f'worker_{worker}{suffix}') != cfg['v2_success_outputs_sha256'][str(worker)][suffix]:
                    raise ValueError(f'Original successful worker {worker} changed')
        workers.append({'path': path, 'sha256': sha(path), 'data': json.loads(path.read_text())})
    # Both combined tables were written before the original aggregation path
    # error. Their bytes must be reproduced by the unchanged aggregator.
    partial = {name: sha(V4_ROOT / 'results' / name) for name in ('pilot_aggregate.tsv', 'pilot_vs_canonical.tsv')}
    OUT.mkdir(exist_ok=False)
    (OUT / 'materialization.receipt.json').symlink_to(V2_ROOT / 'materialization.receipt.json')
    results = OUT / 'results'
    results.mkdir(exist_ok=False)
    for worker in range(1, 5):
        source_root = V4_ROOT if worker == 2 else V2_ROOT
        for suffix in ('.tsv', '.summary.json', '.log'):
            (results / f'worker_{worker}{suffix}').symlink_to(source_root / 'results' / f'worker_{worker}{suffix}')
    mod.PILOT = CONFIG
    mod.RULE = ROOT / cfg['rule_relative']
    mod.CANONICAL = Path(cfg['canonical_aggregate_path'])
    result = mod.aggregate(cfg, rule, OUT, workers)
    for name, digest in partial.items():
        if sha(results / name) != digest:
            raise ValueError(f'Frozen partial {name} differs from repaired aggregation')
    receipt = {'status': result['status'], 'interpretation': 'AGGREGATION_PATH_REPAIR_ONLY',
               'v4_worker_exit_sha256': sha(V4_ROOT / 'results' / 'exit.json'),
               'v4_worker_2_tsv_sha256': exit_info['worker_2_tsv_sha256'],
               'v4_partial_output_sha256': partial,
               'decision_sha256': sha(results / 'pilot_decision.json'),
               'worker_reruns': 0, 'original_v2_rule_sha256': cfg['rule_sha256']}
    mod.save_new(OUT / 'repair_receipt.json', receipt)
    print(json.dumps({'status': result['status'], 'pilot_tested': result['pilot_tested'],
                      'pilot_not_run': result['pilot_status_counts'].get('NOT_RUN', 0),
                      'pilot_failed': result['pilot_status_counts'].get('FAILED', 0),
                      'canonical_tested': result['canonical_tested'],
                      'advance': result['advance_to_full_trait_only_screen']}, sort_keys=True))


if __name__ == '__main__':
    main()

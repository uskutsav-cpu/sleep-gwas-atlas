#!/usr/bin/env python3
"""Receipt-bound v4 repair: reuse v2 successes and rerun only worker 2."""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / 'brain6/config/confirmatory_source_rescue_20260927/insomnia_native_n_pilot_v4.json'
V2_CODE = ROOT / 'brain6/scripts/run_insomnia_native_n_pilot_v2.py'
V4_R = ROOT / 'brain6/scripts/run_insomnia_native_n_pilot_v4.R'
CORE = ROOT / 'brain6/scripts/insomnia_native_n_pilot_v4_core.R'
V2_ROOT = Path('/Volumes/Extreme SSD/brain6-work/confirmatory-source-rescue-20260927/insomnia-native-n-pilot-v2')
V4_ROOT = Path('/Volumes/Extreme SSD/brain6-work/confirmatory-source-rescue-20260927/insomnia-native-n-pilot-v4')
RSCRIPT = Path('/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin/Rscript')


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def save_new(path: Path, value: object) -> None:
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, sort_keys=True, indent=2)
        stream.write('\n')


def check(path: Path, expected: str) -> None:
    observed = sha(path)
    if observed != expected:
        raise ValueError(f'Pinned hash mismatch: {path}: {observed}')


def main() -> None:
    cfg = json.loads(CONFIG.read_text())
    if cfg['analysis_id'] != 'brain6_insomnia_native_n_trait_only_pilot_v4' or cfg['worker_to_run'] != 2:
        raise ValueError('Unexpected repair identity or worker')
    if cfg['v2_root'] != str(V2_ROOT) or cfg['output_root'] != str(V4_ROOT) or V4_ROOT.exists():
        raise FileExistsError('Output root exists or frozen root changed; refuse restart')
    pinned = {
        Path(__file__): cfg['controller_sha256'],
        V2_CODE: cfg['v2_python_runner_sha256'], V4_R: cfg['r_runner_sha256'],
        CORE: cfg['core_sha256'], RSCRIPT: cfg['rscript_sha256'],
        ROOT / cfg['v2_config_relative']: cfg['v2_config_sha256'],
        ROOT / cfg['rule_relative']: cfg['rule_sha256'],
        ROOT / cfg['protocol_relative']: cfg['protocol_sha256'],
        V2_ROOT / 'materialization.receipt.json': cfg['materialization_receipt_sha256'],
        V2_ROOT / 'selected.loci': cfg['selected_loci_sha256'],
        V2_ROOT / 'results' / 'exit.json': cfg['v2_exit_sha256'],
        ROOT / cfg['selection_table_relative']: cfg['selection_table_sha256'],
        ROOT / cfg['selection_receipt_relative']: cfg['selection_receipt_sha256'],
        ROOT / cfg['loci_file_relative']: cfg['loci_file_sha256'],
        Path(cfg['canonical_provenance_path']): cfg['canonical_provenance_sha256'],
        Path(cfg['canonical_aggregate_path']): cfg['canonical_aggregate_sha256'],
        Path(cfg['reference_provenance_path']): cfg['reference_provenance_sha256'],
        Path(cfg['source_readme_path']): cfg['source_readme_sha256'],
        Path(cfg['source_path']): cfg['source_sha256'],
    }
    for path, digest in pinned.items():
        check(path, digest)
    v2 = json.loads((ROOT / cfg['v2_config_relative']).read_text())
    receipt = json.loads((V2_ROOT / 'materialization.receipt.json').read_text())
    rule = json.loads((ROOT / cfg['rule_relative']).read_text())
    if receipt['config_sha256'] != cfg['v2_config_sha256'] or v2['locus_ids'] != cfg['locus_ids']:
        raise ValueError('V2 receipt or locus identity mismatch')
    if rule['advance_if_any'] != cfg['advance_if_any'] or rule['canonical_aggregate_sha256'] != cfg['canonical_aggregate_sha256']:
        raise ValueError('Original advancement gate changed')
    if len(cfg['locus_ids']) != 88 or len(set(cfg['locus_ids'])) != 88 or v2['worker_count'] != 4:
        raise ValueError('Pilot cohort or worker count changed')
    if cfg['source_sha256'] != v2['source_sha256'] or cfg['source_n_semantics'] != v2['source_n_semantics']:
        raise ValueError('Source semantics changed')
    if set(receipt['records']) != set(cfg['locus_ids']):
        raise ValueError('Materialization receipt lacks selected loci')
    workers = []
    found = []
    for worker in range(1, 5):
        path = V2_ROOT / f'worker_{worker}.config.json'
        check(path, cfg['v2_worker_configs_sha256'][str(worker)])
        data = json.loads(path.read_text())
        if data['worker_id'] != worker or data['pilot_config_sha256'] != cfg['v2_config_sha256']:
            raise ValueError('V2 worker configuration changed')
        found.extend(data['locus_ids'])
        workers.append({'path': path, 'sha256': sha(path), 'data': data})
    if len(found) != 88 or len(set(found)) != 88 or set(found) != set(cfg['locus_ids']):
        raise ValueError('V2 worker partitions changed')
    for locus, item in receipt['records'].items():
        directory = V2_ROOT / 'inputs' / f'locus_{locus}'
        check(directory / 'input_info.tsv', item['input_info_sha256'])
        check(directory / 'insomnia.sumstats.tsv.gz', item['sumstats_sha256'])
        if item['rows'] != item['original_rows']:
            raise ValueError('Original vs native-N SNP row count changed')
    for worker in (1, 3, 4):
        for suffix in ('.tsv', '.summary.json', '.log'):
            check(V2_ROOT / 'results' / f'worker_{worker}{suffix}', cfg['v2_success_outputs_sha256'][str(worker)][suffix])
        summary = json.loads((V2_ROOT / 'results' / f'worker_{worker}.summary.json').read_text())
        if summary['rows'] != 22 or summary['failed'] != 0:
            raise ValueError(f'V2 worker {worker} not complete')
    # All checks finish before creating the exclusive v4 output directory.
    V4_ROOT.mkdir(exist_ok=False)
    results = V4_ROOT / 'results'
    results.mkdir(exist_ok=False)
    for worker in (1, 3, 4):
        for suffix in ('.tsv', '.summary.json', '.log'):
            (results / f'worker_{worker}{suffix}').symlink_to(V2_ROOT / 'results' / f'worker_{worker}{suffix}')
    worker2 = workers[1]['data'].copy()
    worker2.update({
        'analysis_id': cfg['analysis_id'], 'scope': cfg['scope'],
        'core_path': str(CORE),
        'input_rows': {loc: receipt['records'][loc]['rows'] for loc in workers[1]['data']['locus_ids']},
        'output_tsv': str(results / 'worker_2.tsv'),
        'summary_json': str(results / 'worker_2.summary.json'),
    })
    worker2_config = V4_ROOT / 'worker_2.config.json'
    save_new(worker2_config, worker2)
    log_path = results / 'worker_2.log'
    env = dict(os.environ)
    env.update({'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1', 'MKL_NUM_THREADS': '1'})
    command = [str(RSCRIPT), str(V4_R), str(worker2_config)]
    launch = {'analysis_id': cfg['analysis_id'], 'v4_config_sha256': sha(CONFIG),
              'v2_success_workers_reused': [1, 3, 4], 'worker_rerun': 2,
              'worker_config_sha256': sha(worker2_config), 'command': command}
    save_new(results / 'launch.json', launch)
    print(json.dumps({'status': 'WORKER_2_RUNNING', 'reused_workers': [1, 3, 4], 'output_root': str(V4_ROOT)}), flush=True)
    with log_path.open('xb') as log:
        code = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT, check=False).returncode
    exit_info = {'analysis_id': cfg['analysis_id'], 'worker_2_exit_code': code,
                 'worker_2_log_sha256': sha(log_path), 'launch_sha256': sha(results / 'launch.json')}
    if (results / 'worker_2.tsv').exists():
        exit_info['worker_2_tsv_sha256'] = sha(results / 'worker_2.tsv')
    if (results / 'worker_2.summary.json').exists():
        exit_info['worker_2_summary_sha256'] = sha(results / 'worker_2.summary.json')
    save_new(results / 'exit.json', exit_info)
    if code:
        raise RuntimeError(f'Worker 2 still failed; inspect receipt: {code}')
    spec = importlib.util.spec_from_file_location('insomnia_v2_runner', V2_CODE)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.PILOT = CONFIG
    module.RULE = ROOT / cfg['rule_relative']
    module.CANONICAL = Path(cfg['canonical_aggregate_path'])
    # The aggregation code is unchanged. It checks all 88 statuses and the
    # original gate; v4 provenance explicitly points to reused v2 workers.
    outcome = module.aggregate(cfg, rule, V4_ROOT, workers)
    outcome['v2_success_workers_reused'] = [1, 3, 4]
    outcome['worker_2_rerun_only'] = True
    save_new(results / 'v4_receipt.json', outcome)
    print(json.dumps({'status': outcome['status'], 'pilot_tested': outcome['pilot_tested'],
                      'canonical_tested': outcome['canonical_tested'],
                      'advance': outcome['advance_to_full_trait_only_screen']}), flush=True)


if __name__ == '__main__':
    main()

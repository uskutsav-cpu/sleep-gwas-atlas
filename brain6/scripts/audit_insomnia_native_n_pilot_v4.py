#!/usr/bin/env python3
"""Audit the fixed 88-locus pilot and its non-promotional family projection."""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'brain6/results/confirmatory_source_rescue_20260927/other_traits'
V2 = Path('/Volumes/Extreme SSD/brain6-work/confirmatory-source-rescue-20260927/insomnia-native-n-pilot-v2')
V4 = Path('/Volumes/Extreme SSD/brain6-work/confirmatory-source-rescue-20260927/insomnia-native-n-pilot-v4')
REPAIR = V4 / 'aggregation_repair_v1'
CANON = ROOT / 'work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/results/canonical_family_results.tsv'
LOCI = ROOT / 'ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile'
CAND = ROOT / 'brain6/results/loci/placo_five_track_candidate_loci_v1/candidate_loci.tsv'
OUT = BASE / 'insomnia_native_n_pilot_v4_audit'


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def table(path: Path, delimiter: str = '\t') -> list[dict[str, str]]:
    with path.open(newline='') as f:
        return list(csv.DictReader(f, delimiter=delimiter, skipinitialspace=True))


def main() -> None:
    if OUT.exists():
        raise FileExistsError('Audit already exists; refuse overwrite')
    cfg_path = ROOT / 'brain6/config/confirmatory_source_rescue_20260927/insomnia_native_n_pilot_v4.json'
    cfg = json.loads(cfg_path.read_text())
    if sha(cfg_path) != '9cd1f0a18a17e370e68e1e6170f8c9edda1a7b7ad86e5887505bd46e0d14d47c':
        raise ValueError('Frozen v4 config changed')
    rule_path = ROOT / cfg['rule_relative']
    if sha(rule_path) != cfg['rule_sha256']:
        raise ValueError('Frozen advancement gate changed')
    decision_path = REPAIR / 'results' / 'pilot_decision.json'
    decision = json.loads(decision_path.read_text())
    repair_receipt = json.loads((REPAIR / 'repair_receipt.json').read_text())
    if repair_receipt['decision_sha256'] != sha(decision_path) or repair_receipt['worker_reruns'] != 0:
        raise ValueError('Aggregation repair receipt mismatch')
    exit_info = json.loads((V4 / 'results' / 'exit.json').read_text())
    if exit_info['worker_2_exit_code'] != 0 or sha(V4 / 'results' / 'worker_2.tsv') != exit_info['worker_2_tsv_sha256']:
        raise ValueError('Worker 2 completion receipt mismatch')
    rows = table(REPAIR / 'results' / 'pilot_vs_canonical.tsv')
    if len(rows) != 88 or {r['locus_id'] for r in rows} != set(cfg['locus_ids']):
        raise ValueError('Pilot family changed')
    transitions = Counter((r['canonical_status'], r['pilot_status']) for r in rows)
    if transitions != {('TESTED', 'TESTED'): 44, ('NOT_RUN', 'NOT_RUN'): 41, ('NOT_RUN', 'TESTED'): 3}:
        raise ValueError('Unexpected pilot transitions')
    if decision['canonical_tested'] != 44 or decision['pilot_tested'] != 47 or decision['pilot_status_counts'] != {'TESTED': 47, 'NOT_RUN': 41}:
        raise ValueError('Pilot decision/count mismatch')
    if decision['advance_to_full_trait_only_screen'] or decision['family_promotion_authorized']:
        raise ValueError('Advancement state changed')
    loci = table(LOCI, ' ')
    by_id = {r['LOC']: r for r in loci}
    candidates = table(CAND)
    excluded = {l['LOC'] for l in loci for c in candidates if l['CHR'] == c['CHR']
                and int(l['START']) <= int(c['STOP']) and int(c['START']) <= int(l['STOP'])}
    if len(excluded) != 38:
        raise ValueError('Candidate-excluded pilot design changed')
    eligible = Counter()
    for row in table(CANON):
        if row['phen'] == 'insomnia' and row['status'] == 'NOT_RUN' and row['locus_id'] not in excluded:
            eligible[by_id[row['locus_id']]['CHR']] += 1
    converted = Counter(by_id[r['locus_id']]['CHR'] for r in rows
                        if r['canonical_status'] == 'NOT_RUN' and r['pilot_status'] == 'TESTED')
    if sum(eligible.values()) != 670 or sum(converted.values()) != 3:
        raise ValueError('Projection strata changed')
    diagnostic_estimate = sum(eligible[ch] * converted[ch] / 2 for ch in eligible)
    # This estimates candidate-excluded insomnia cells under the two-per-chromosome
    # hash sample. It is not a validated full-screen or family recovery count.
    summary = {
        'analysis_id': cfg['analysis_id'], 'status': decision['status'],
        'canonical': {'tested': 44, 'not_run': 44, 'failed': 0},
        'pilot': {'tested': 47, 'not_run': 41, 'failed': 0},
        'tested_delta': 3, 'absolute_tested_fraction_gain': decision['absolute_tested_fraction_gain'],
        'canonical_low_h2_not_run': decision['canonical_low_local_h2_not_run'],
        'pilot_low_h2_not_run': decision['pilot_low_local_h2_not_run'],
        'relative_low_h2_reduction': decision['relative_low_local_h2_not_run_reduction'],
        'conversion_locus_ids': [r['locus_id'] for r in rows if r['canonical_status'] == 'NOT_RUN' and r['pilot_status'] == 'TESTED'],
        'empty_locus_id': '950', 'empty_locus_status': 'NOT_RUN_NO_SNPS',
        'candidate_excluded_not_run_cells': sum(eligible.values()),
        'candidate_excluded_diagnostic_recovery_estimate': diagnostic_estimate,
        'diagnostic_projection_method': 'Within each chromosome, candidate-excluded canonical NOT_RUN count times observed conversions among two hash-ranked sampled NOT_RUN loci divided by two; not a full-family estimator or decision metric.',
        'validated_family_recovery_cells': 0, 'canonical_family_not_run': 3720,
        'validated_projected_family_not_run': 3720, 'frozen_family_ceiling': 873,
        'full_insomnia_screen_authorized': False, 'confirmatory_source_admitted': False,
        'source_status': 'ADMISSIBLE_SENSITIVITY_ONLY',
        'source_sha256': {str(p): sha(p) for p in (cfg_path, rule_path, decision_path,
                          REPAIR / 'repair_receipt.json', V4 / 'results' / 'exit.json',
                          V4 / 'results' / 'worker_2.tsv', V2 / 'materialization.receipt.json',
                          CANON, LOCI, CAND)},
        'builder_sha256': sha(Path(__file__)),
    }
    OUT.mkdir(parents=True)
    with (OUT / 'summary.json').open('x') as f:
        json.dump(summary, f, indent=2, sort_keys=True)
        f.write('\n')
    print(json.dumps({'status': summary['status'], 'validated_family_recovery': 0,
                      'diagnostic_candidate_excluded_recovery_estimate': diagnostic_estimate}))


if __name__ == '__main__':
    main()

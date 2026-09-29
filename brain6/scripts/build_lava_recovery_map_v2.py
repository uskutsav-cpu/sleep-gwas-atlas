#!/usr/bin/env python3
"""Read-only locus and candidate overlap map for the frozen LAVA family."""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CANON = ROOT / 'work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/results/canonical_family_results.tsv'
LOCI = ROOT / 'ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile'
CAND = ROOT / 'brain6/results/loci/placo_five_track_candidate_loci_v1/candidate_loci.tsv'
OUT = ROOT / 'brain6/results/confirmatory_source_rescue_20260927/phase2_recovery_map_v2'


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open('x', newline='') as f:
        w = csv.DictWriter(f, delimiter='\t', fieldnames=fields, lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    if OUT.exists():
        raise FileExistsError('Refuse to overwrite recovery map')
    with LOCI.open() as f:
        rd = csv.DictReader(f, delimiter=' ', skipinitialspace=True)
        blocks = [{k: int(v) for k, v in row.items()} for row in rd]
    if len(blocks) != 2495 or {r['LOC'] for r in blocks} != set(range(1, 2496)):
        raise ValueError('LAVA block universe changed')
    canon = read_tsv(CANON)
    if len(canon) != 17465 or Counter(r['status'] for r in canon) != {'TESTED': 13745, 'NOT_RUN': 3720}:
        raise ValueError('Canonical LAVA family changed')
    index = {(r['phen'], int(r['locus_id'])): r for r in canon}
    if len(index) != len(canon):
        raise ValueError('Duplicate canonical trait/locus')
    candidates = read_tsv(CAND)
    if len(candidates) != 25 or len({r['locus_id'] for r in candidates}) != 25:
        raise ValueError('Protected candidate table changed')
    overlaps = []
    summaries = []
    for candidate in candidates:
        pair = candidate['pair_id']
        left, right = pair.split('__')
        chrom, start, stop = (int(candidate[k]) for k in ('CHR', 'START', 'STOP'))
        hit = [b for b in blocks if b['CHR'] == chrom and b['START'] <= stop and b['STOP'] >= start]
        if not hit:
            raise ValueError(f'Candidate has no LAVA block: {candidate["locus_id"]}')
        tested_both = 0
        for b in hit:
            a, z = index[(left, b['LOC'])], index[(right, b['LOC'])]
            tested_both += a['status'] == z['status'] == 'TESTED'
            overlaps.append({'candidate_locus_id': candidate['locus_id'], 'pair_id': pair,
                             'candidate_chr': chrom, 'candidate_start': start, 'candidate_stop': stop,
                             'lead_variant': candidate['lead_variant'], 'lava_locus_id': b['LOC'],
                             'lava_start': b['START'], 'lava_stop': b['STOP'],
                             'left_trait': left, 'left_status': a['status'], 'left_reason': a['reason'],
                             'left_p': a['p'], 'right_trait': right, 'right_status': z['status'],
                             'right_reason': z['reason'], 'right_p': z['p'],
                             'both_univariate_tested': a['status'] == z['status'] == 'TESTED'})
        summaries.append({'candidate_locus_id': candidate['locus_id'], 'pair_id': pair,
                          'lead_variant': candidate['lead_variant'], 'overlapping_lava_loci': len(hit),
                          'both_univariate_tested_loci': tested_both,
                          'blocked_univariate_loci': len(hit) - tested_both,
                          'interpretation': 'READ_ONLY_CANONICAL_UNIVARIATE;PAIRWISE_NOT_EXECUTED'})
    OUT.mkdir(parents=True)
    write_tsv(OUT / 'candidate_lava_locus_overlap.tsv', list(overlaps[0]), overlaps)
    write_tsv(OUT / 'candidate_recovery_summary.tsv', list(summaries[0]), summaries)
    receipt = {'schema_version': 1, 'scope': 'READ_ONLY_FROZEN_CANONICAL_RECOVERY_MAP',
               'source_sha256': {str(p.relative_to(ROOT)): sha(p) for p in (CANON, LOCI, CAND)},
               'builder_sha256': sha(Path(__file__)), 'candidate_count': len(summaries),
               'overlap_rows': len(overlaps), 'pair_count': len({r['pair_id'] for r in summaries}),
               'outputs_sha256': {p.name: sha(p) for p in OUT.glob('*.tsv')}}
    with (OUT / 'receipt.json').open('x') as f:
        json.dump(receipt, f, indent=2, sort_keys=True)
        f.write('\n')
    print(json.dumps({'candidate_count': len(summaries), 'overlap_rows': len(overlaps),
                      'blocked_candidates': sum(r['blocked_univariate_loci'] > 0 for r in summaries)}))


if __name__ == '__main__':
    main()

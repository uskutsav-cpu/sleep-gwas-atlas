#!/usr/bin/env python3
"""Project the frozen LAVA family after all source-first decisions to date."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'brain6/results/confirmatory_source_rescue_20260927'
TRAITS = BASE / 'phase1/trait_failure_burden.tsv'
INSOMNIA = BASE / 'other_traits/insomnia_native_n_pilot_v4_audit/summary.json'
LONGSLEEP = BASE / 'public_longsleep/SOURCE_EVIDENCE_LEDGER.md'
OTHER = BASE / 'other_traits/SOURCE_FIRST_LEDGER.md'
OUT = BASE / 'phase5_family_projection_v2'
STATUS = {
    'longsleep': ('EXACT_SOURCE_METADATA_AND_MODEL_HOLD', 'No source-admitted exact N/model mapping; prior pilot did not advance'),
    'insomnia': ('NATIVE_N_PILOT_GATE_FAILED', 'V4 technical pilot improved 3/88 but missed both original advancement gates'),
    'parkinson': ('METHOD_HOLD', 'Mixed direct/proxy cases and variant count semantics unresolved'),
    'mdd': ('SOURCE_N_MISSING_OR_REPLACEMENT_HOLD', 'Original lacks N; newer source changes cohort/model and does not pass admission'),
    'adhd': ('METHOD_HOLD', 'Varying per-variant case fractions and cohort mix unresolved'),
    'bipolar': ('METHOD_HOLD', 'Released count fields are effective, not verified literal analyzed counts'),
    'scz': ('METHOD_HOLD', 'Effective-N convention and per-variant literal counts unresolved'),
}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    if OUT.exists():
        raise FileExistsError('Projection exists; refuse overwrite')
    with TRAITS.open(newline='') as f:
        traits = list(csv.DictReader(f, delimiter='\t'))
    if len(traits) != 7 or sum(int(t['NOT_RUN']) for t in traits) != 3720:
        raise ValueError('Canonical trait failure burden changed')
    insomnia = json.loads(INSOMNIA.read_text())
    if insomnia['status'] != 'NO_MATERIAL_PILOT_IMPROVEMENT' or insomnia['validated_family_recovery_cells'] != 0:
        raise ValueError('Insomnia terminal admission changed')
    rows = []
    cumulative_perfect = 0
    for trait in traits:
        name = trait['trait_id']
        burden = int(trait['NOT_RUN'])
        cumulative_perfect += burden
        rows.append({'rank': trait['rank_by_not_run'], 'trait_id': name,
                     'source_gwas': trait['implicated_source'],
                     'canonical_tested': trait['TESTED'], 'canonical_not_run': burden,
                     'maximum_theoretical_recovery_cells': burden,
                     'source_gate': STATUS[name][0], 'source_gate_reason': STATUS[name][1],
                     'validated_recovery_cells': 0,
                     'validated_projected_family_not_run': 3720,
                     'hypothetical_residual_if_all_ranked_traits_to_here_perfectly_repaired': 3720 - cumulative_perfect,
                     'frozen_family_ceiling': 873})
    if [r['trait_id'] for r in rows[:4]] != ['longsleep', 'insomnia', 'parkinson', 'mdd'] or rows[2]['hypothetical_residual_if_all_ranked_traits_to_here_perfectly_repaired'] != 1122 or rows[3]['hypothetical_residual_if_all_ranked_traits_to_here_perfectly_repaired'] != 533:
        raise ValueError('Frozen minimal perfect-repair arithmetic changed')
    OUT.mkdir(parents=True)
    table = OUT / 'family_projection.tsv'
    with table.open('x', newline='') as f:
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    summary = {'schema_version': 1, 'canonical_tested': 13745, 'canonical_not_run': 3720,
               'canonical_failed': 0, 'frozen_maximum_not_run': 873,
               'minimum_required_validated_recovery': 2847,
               'validated_recovery_to_date': 0, 'validated_projected_not_run': 3720,
               'smallest_hypothetical_perfect_repair_set': ['longsleep', 'insomnia', 'parkinson', 'mdd'],
               'perfect_repair_residual_after_three': 1122, 'perfect_repair_residual_after_four': 533,
               'full_rescue_launch_authorized': False,
               'source_sha256': {str(p.relative_to(ROOT)): sha(p) for p in (TRAITS, INSOMNIA, LONGSLEEP, OTHER)},
               'projection_tsv_sha256': sha(table), 'builder_sha256': sha(Path(__file__))}
    with (OUT / 'summary.json').open('x') as f:
        json.dump(summary, f, indent=2, sort_keys=True)
        f.write('\n')
    print(json.dumps({'validated_projected_not_run': 3720, 'ceiling': 873,
                      'smallest_perfect_repair_set_size': 4}))


if __name__ == '__main__':
    main()

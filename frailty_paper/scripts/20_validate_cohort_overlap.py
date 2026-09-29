#!/usr/bin/env python3
"""Validate the Phase 8 cohort ledger against the frozen analysis panel."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path

SLEEP = {
    'insomnia', 'sleepdur', 'shortsleep', 'longsleep', 'chronotype', 'sleepiness',
    'napping', 'snoring', 'sleep_apnea', 'sleep_efficiency',
    'accel_sleep_duration', 'sleep_timing',
}
REQUIRED = {
    'trait_id', 'role', 'source_id', 'accession', 'source_cohort', 'uk_biobank',
    'finngen', '23andme', 'charge', 'consortium_cohorts', 'other_known_cohorts',
    'source_status', 'exact_participant_intersection', 'evidence_file', 'notes',
}
ALLOWED = {'YES', 'NO', 'UNKNOWN', 'NOT_APPLICABLE'}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def validate(repo: Path) -> list[str]:
    ledger = read_tsv(repo / 'frailty_analysis/manifests/cohort_overlap.tsv')
    panel = {r['trait_id']: r['source_id'] for r in read_tsv(repo / 'config/analysis_panel.tsv') if r['trait_id'] in SLEEP}
    issues: list[str] = []
    if len(ledger) < len(SLEEP):
        issues.append(f'expected at least {len(SLEEP)} registered endpoint rows; found {len(ledger)}')
    if not ledger or not REQUIRED.issubset(ledger[0]):
        issues.append(f'missing required columns: {sorted(REQUIRED - set(ledger[0] if ledger else []))}')
        return issues
    sleep_rows = [r for r in ledger if r['role'] == 'primary_sleep']
    observed = {r['trait_id']: r for r in sleep_rows}
    if set(observed) != SLEEP:
        issues.append(f'sleep panel mismatch; missing={sorted(SLEEP-set(observed))}, extra={sorted(set(observed)-SLEEP)}')
    for trait in SLEEP & set(observed):
        if observed[trait]['source_id'] != panel.get(trait):
            issues.append(f'{trait}: source_id does not match frozen analysis_panel.tsv')
    for i, row in enumerate(ledger, start=2):
        for key in ('uk_biobank', 'finngen', '23andme', 'charge'):
            if row[key] not in ALLOWED:
                issues.append(f'row {i}: invalid {key}={row[key]!r}')
        if row['exact_participant_intersection'] != 'UNKNOWN':
            issues.append(f'row {i}: exact participant intersection must remain UNKNOWN absent verified counts')
        if not row['evidence_file'] or not row['notes']:
            issues.append(f'row {i}: evidence path and limitations note are required')
    if len(observed) != len(sleep_rows):
        issues.append('duplicate sleep trait rows')
    return issues


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    issues = validate(args.repo)
    if issues:
        for issue in issues:
            print(f'ERROR: {issue}')
        return 1
    print(f'cohort ledger valid: {len(read_tsv(args.repo / "frailty_analysis/manifests/cohort_overlap.tsv"))} rows; frozen 12-trait panel/source IDs match; exact intersections unknown')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())

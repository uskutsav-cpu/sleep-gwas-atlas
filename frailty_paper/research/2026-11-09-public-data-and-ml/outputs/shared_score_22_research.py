#!/usr/bin/env python3
"""Compute an unvalidated 22-marker shared-direction allele count, research only."""
import argparse
import csv
import sys
from pathlib import Path
from score_22_research import calculate, read_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dosages', required=True, help='Long-format CSV documented in the score README')
    parser.add_argument('--output', required=True, help='New output file; never overwrites')
    parser.add_argument('--allow-verified-reversal', action='store_true', help='Permit 2 minus dosage only after independent verification of every reversed biallelic diploid contrast')
    args = parser.parse_args()
    try:
        model = read_rows(Path(__file__).with_name('unvalidated_22_shared_direction_score.csv'))
        if len(model) != 22 or any(row['weight'] != '1' for row in model):
            raise ValueError('Expected 22 unit-weight variants')
        compatible = [{**row, 'FI_beta_weight': 1, 'FFS_beta_weight_inferred': 1,
                       'UKB_insomnia_beta_weight': 1} for row in model]
        result = calculate(compatible, read_rows(args.dosages), args.allow_verified_reversal)
        with Path(args.output).open('x', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=['sample_id', 'n_variants', 'shared_direction_allele_count'])
            writer.writeheader()
            for row in result:
                writer.writerow({'sample_id': row['sample_id'], 'n_variants': 22,
                                 'shared_direction_allele_count': row['FI_raw_sum']})
    except (KeyError, ValueError, OSError) as exc:
        print(f'Scoring failed: {exc}', file=sys.stderr)
        return 1
    print(f'Wrote unvalidated research counts for {len(result)} samples.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

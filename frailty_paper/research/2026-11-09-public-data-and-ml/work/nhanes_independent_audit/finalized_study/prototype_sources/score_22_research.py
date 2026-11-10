#!/usr/bin/env python3
"""Compute unvalidated research-only 22-marker FI, FFS, and insomnia sums.

Dosages are preprocessed and must represent the listed allele pair on a
verified biallelic diploid contrast. See score_22_research_readme.md.
"""
import argparse
import csv
import math
import sys
from pathlib import Path

WEIGHTS = ('FI_beta_weight', 'FFS_beta_weight_inferred', 'UKB_insomnia_beta_weight')


def read_rows(path):
    with Path(path).open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        if not fields or any(not name or not name.strip() for name in fields) or len(fields) != len(set(fields)):
            raise ValueError('CSV requires nonempty, unique column names')
        rows = list(reader)
        if any(None in row or any(value is None for value in row.values()) for row in rows):
            raise ValueError('CSV row has too many or too few fields')
        return rows


def calculate(weights, observations, allow_verified_reversal=False):
    model = {}
    for row in weights:
        rs = row['variant']
        a, b = row['counted_allele'].upper(), row['other_allele'].upper()
        if not rs or rs in model or a not in ('A', 'C', 'G', 'T') or b not in ('A', 'C', 'G', 'T') or a == b:
            raise ValueError(f'Invalid or duplicate model variant: {rs}')
        coefficients = [float(row[k]) for k in WEIGHTS]
        if not all(math.isfinite(x) for x in coefficients):
            raise ValueError(f'Invalid weight: {rs}')
        model[rs] = (a, b, coefficients)
    if len(model) != 22:
        raise ValueError('Expected exactly 22 unique prototype variants')
    grouped = {}
    for row in observations:
        sample, rs = row['sample_id'], row['variant']
        if not sample or rs not in model:
            raise ValueError(f'Empty sample ID or variant not in model: {rs}')
        variants = grouped.setdefault(sample, {})
        if rs in variants:
            raise ValueError(f'Duplicate observation for sample {sample}, variant {rs}')
        a, b, _ = model[rs]
        da, db = row['dosage_allele'].upper(), row['other_allele'].upper()
        try:
            dosage = float(row['dosage'])
        except ValueError:
            raise ValueError(f'Missing or invalid dosage for sample {sample}, variant {rs}') from None
        if not math.isfinite(dosage) or not 0 <= dosage <= 2:
            raise ValueError(f'Dosage outside [0, 2] for sample {sample}, variant {rs}')
        if (da, db) == (a, b):
            aligned = dosage
        elif (da, db) == (b, a):
            if not allow_verified_reversal:
                raise ValueError(f'Reversed dosage for sample {sample}, variant {rs}; provide counted-allele dosage or explicitly verify biallelic diploid reversal')
            aligned = 2 - dosage
        else:
            raise ValueError(f'Allele mismatch for sample {sample}, variant {rs}')
        variants[rs] = aligned
    if not grouped:
        raise ValueError('No observations provided')
    result = []
    for sample in sorted(grouped):
        variants = grouped[sample]
        missing = set(model) - variants.keys()
        if missing:
            raise ValueError(f'Sample {sample} is missing {len(missing)} model variants')
        sums = [math.fsum(variants[rs] * model[rs][2][i] for rs in model) for i in range(3)]
        result.append({'sample_id': sample, 'n_variants': 22, 'FI_raw_sum': sums[0], 'FFS_raw_sum': sums[1], 'UKB_insomnia_raw_sum': sums[2]})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dosages', required=True, help='Long-format CSV, see README')
    parser.add_argument('--weights', default=str(Path(__file__).with_name('unvalidated_22_variant_score_prototype.csv')))
    parser.add_argument('--output', required=True, help='New file path; an existing output is never overwritten')
    parser.add_argument('--allow-verified-reversal', action='store_true', help='Permit 2 minus dosage only after independent verification of every reversed biallelic diploid contrast')
    args = parser.parse_args()
    try:
        results = calculate(read_rows(args.weights), read_rows(args.dosages), args.allow_verified_reversal)
        with Path(args.output).open('x', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
    except (KeyError, ValueError, OSError) as exc:
        print(f'Scoring failed: {exc}', file=sys.stderr)
        return 1
    print(f'Wrote unvalidated research sums for {len(results)} samples.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

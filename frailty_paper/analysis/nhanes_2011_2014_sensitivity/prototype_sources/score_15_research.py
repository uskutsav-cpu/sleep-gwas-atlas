#!/usr/bin/env python3
"""Calculate three unvalidated 15-marker sums from QCed research dosages.

This does not predict individual risk. Inputs must already be verified as
biallelic, diploid and on the model's genome build and strand.
"""
import argparse
import csv
import math
import sys
from pathlib import Path
from score_22_research import read_rows

WEIGHTS = ('FI_beta_weight', 'FFS_beta_weight_inferred', 'insomnia_beta_weight')

def calculate(weights, observations, allow_verified_reversal=False):
    model = {}
    for row in weights:
        rs = row['variant']
        a, b = row['counted_allele'].upper(), row['other_allele'].upper()
        if not rs or rs in model or a not in 'ACGT' or b not in 'ACGT' or len(a) != 1 or len(b) != 1 or a == b:
            raise ValueError(f'Invalid or duplicate model variant: {rs}')
        v = tuple(float(row[k]) for k in WEIGHTS)
        if not all(math.isfinite(x) for x in v):
            raise ValueError(f'Invalid weight: {rs}')
        model[rs] = (a, b, v)
    if len(model) != 15:
        raise ValueError('Expected exactly 15 unique prototype variants')
    grouped = {}
    for row in observations:
        sample, rs = row['sample_id'], row['variant']
        if not sample or rs not in model:
            raise ValueError(f'Empty sample ID or variant not in model: {rs}')
        vals = grouped.setdefault(sample, {})
        if rs in vals:
            raise ValueError(f'Duplicate observation for sample {sample}, variant {rs}')
        a, b, _ = model[rs]
        da, db = row['dosage_allele'].upper(), row['other_allele'].upper()
        try:
            dosage = float(row['dosage'])
        except (ValueError, TypeError):
            raise ValueError(f'Missing or invalid dosage for sample {sample}, variant {rs}') from None
        if not math.isfinite(dosage) or not 0 <= dosage <= 2:
            raise ValueError(f'Dosage outside [0, 2] for sample {sample}, variant {rs}')
        if (da, db) == (a, b):
            vals[rs] = dosage
        elif (da, db) == (b, a):
            if not allow_verified_reversal:
                raise ValueError(f'Reversed dosage for sample {sample}, variant {rs}; provide counted-allele dosage or explicitly verify biallelic diploid reversal')
            vals[rs] = 2 - dosage
        else:
            raise ValueError(f'Allele mismatch for sample {sample}, variant {rs}')
    if not grouped:
        raise ValueError('No observations provided')
    result = []
    for sample in sorted(grouped):
        vals = grouped[sample]
        if set(vals) != set(model):
            raise ValueError(f'Sample {sample} is missing {len(set(model) - set(vals))} model variants')
        sums = [math.fsum(vals[rs] * model[rs][2][i] for rs in model) for i in range(3)]
        result.append({'sample_id': sample, 'n_variants': 15, 'FI_raw_sum': sums[0],
                       'FFS_raw_sum_inferred': sums[1], 'insomnia_raw_sum': sums[2]})
    return result

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dosages', required=True)
    parser.add_argument('--output', required=True, help='New file; existing files are never overwritten')
    parser.add_argument('--weights', default=str(Path(__file__).with_name('unvalidated_15_variant_score_prototype.csv')))
    parser.add_argument('--allow-verified-reversal', action='store_true', help='Permit 2 minus dosage only after independent verification of every reversed biallelic diploid contrast')
    args = parser.parse_args()
    try:
        model = read_rows(args.weights)
        result = calculate(model, read_rows(args.dosages), args.allow_verified_reversal)
        with Path(args.output).open('x', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(result[0]))
            writer.writeheader()
            writer.writerows(result)
    except (KeyError, TypeError, ValueError, OSError) as exc:
        print(f'Scoring failed: {exc}', file=sys.stderr)
        return 1
    print(f'Wrote unvalidated research sums for {len(result)} samples.')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())

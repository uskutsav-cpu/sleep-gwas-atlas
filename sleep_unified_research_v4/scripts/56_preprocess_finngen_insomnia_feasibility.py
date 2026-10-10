#!/usr/bin/env python3
"""Streaming, qualified R13 source diagnostic; no fit or network request."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import math
from pathlib import Path
import re
import sys

from canonical_calibration_common_v4_3 import check_hashes, sha, utc, validate_inherited_lock, write_new

HEADER = ['#chrom', 'pos', 'ref', 'alt', 'rsids', 'nearest_genes', 'pval',
          'mlogp', 'beta', 'sebeta', 'af_alt', 'af_alt_cases', 'af_alt_controls']
VALID = {'A', 'C', 'G', 'T'}
AMBIGUOUS = {frozenset(('A', 'T')), frozenset(('C', 'G'))}
COMP = str.maketrans('ACGT', 'TGCA')
RSID = re.compile(r'(?<![A-Za-z0-9_])rs[0-9]+(?![A-Za-z0-9_])')


def finite(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (ValueError, TypeError, OverflowError):
        return None


def orientation(effect, other, target_a1, target_a2):
    effect, other = effect.upper(), other.upper()
    if effect not in VALID or other not in VALID or effect == other or frozenset((effect, other)) in AMBIGUOUS:
        return None
    if (effect, other) == (target_a1, target_a2):
        return 1
    if (effect, other) == (target_a2, target_a1):
        return -1
    complemented = effect.translate(COMP), other.translate(COMP)
    if complemented == (target_a1, target_a2):
        return 1
    if complemented == (target_a2, target_a1):
        return -1
    return None


def unique_hm3_rsid(value, hm3):
    candidates = set(RSID.findall(value)).intersection(hm3)
    if not candidates:
        return 'not_nonambiguous_hm3', None
    if len(candidates) != 1:
        return 'multiple_hm3_rsids', None
    return None, candidates.pop()


def qualify_row(row, hm3, chain, seen):
    reason, rsid = unique_hm3_rsid(row['rsids'], hm3)
    if reason:
        return reason, None
    if rsid in seen:
        return 'duplicate_retained_hm3_rsid', None
    a1, a2, chromosome, position = hm3[rsid]
    status, target = chain.map_point(row['#chrom'], row['pos'])
    if status != 'mapped':
        return 'liftover_' + status, None
    if target[:2] != (chromosome, position):
        return 'mapped_coordinate_rsid_reference_mismatch', None
    if chromosome == 6 and 25_000_000 <= position <= 34_000_000:
        return 'extended_MHC_GRCh37', None
    effect, other = row['alt'].upper(), row['ref'].upper()
    if target[2] == '-':
        effect, other = effect.translate(COMP), other.translate(COMP)
    sign = orientation(effect, other, a1, a2)
    if sign is None:
        return 'allele_mismatch_or_ambiguous', None
    frequency = finite(row['af_alt'])
    if frequency is None or not 0 <= frequency <= 1 or min(frequency, 1-frequency) <= 0.01:
        return 'maf_at_or_below_0_01_or_invalid', None
    beta, se = finite(row['beta']), finite(row['sebeta'])
    if beta is None or se is None or se <= 0:
        return 'invalid_beta_or_supplied_se', None
    z = beta / se
    if not math.isfinite(z):
        return 'nonfinite_beta_over_se', None
    return 'retained', dict(SNP=rsid, A1=a1, A2=a2, Z=sign*z,
                            source_z=z, supplied_p=finite(row['pval']), reverse_chain=target[2] == '-')


def load_hm3(reference, allele_list):
    alleles = {}
    with Path(allele_list).open(newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            snp, a1, a2 = row['SNP'], row['A1'].upper(), row['A2'].upper()
            if snp in alleles:
                raise RuntimeError('DUPLICATE_REFERENCE_ALLELE_SNP')
            if a1 in VALID and a2 in VALID and a1 != a2 and frozenset((a1, a2)) not in AMBIGUOUS:
                alleles[snp] = a1, a2
    hm3 = {}
    with gzip.open(reference, 'rt', newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            snp = row['SNP']
            if snp not in alleles:
                continue
            if snp in hm3:
                raise RuntimeError('DUPLICATE_REFERENCE_COORDINATE_SNP')
            chromosome, position = int(row['CHR']), int(row['BP'])
            if not 1 <= chromosome <= 22 or position <= 0:
                raise RuntimeError('INVALID_GRCH37_REFERENCE_COORDINATE')
            hm3[snp] = (*alleles[snp], chromosome, position)
    if len(hm3) < 900_000:
        raise RuntimeError('HISTORICAL_NONAMBIGUOUS_HM3_REFERENCE_TOO_SMALL')
    return hm3


def body_hashes(path):
    hashes = {'sha256': hashlib.sha256(), 'md5': hashlib.md5()}
    size = 0
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            size += len(block)
            for digest in hashes.values():
                digest.update(block)
    return dict(bytes=size, **{k:v.hexdigest() for k,v in hashes.items()})


def preprocess(plan_path, plan_hash, lock_fd):
    validate_inherited_lock(lock_fd)
    if sha(plan_path) != plan_hash:
        raise RuntimeError('FROZEN_PLAN_CHANGED')
    plan = json.loads(Path(plan_path).read_text())
    if plan['scope'] != 'FINNGEN_INSOMNIA_SOURCE_FEASIBILITY_ONLY' or plan['new_rg_commands'] != 0:
        raise RuntimeError('SOURCE_FEASIBILITY_SCOPE_CHANGED')
    check_hashes(plan['dependencies_sha256'])
    source, output, receipt = map(Path, [plan['source'], plan['derivative'], plan['preprocessing_receipt']])
    expected = plan['source_identity']
    before = body_hashes(source)
    if before != expected:
        raise RuntimeError('RAW_BODY_IDENTITY_CHANGED')
    if output.exists() or receipt.exists():
        raise RuntimeError('PRIOR_FEASIBILITY_ATTEMPT_PRESERVED_NO_OVERWRITE')
    spec = importlib.util.spec_from_file_location('frozen_finngen_liftover', plan['liftover_code'])
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    chain, chain_metadata = module.load_chain(plan['chain'], expected_sha256=plan['chain_sha256'], expected_bytes=plan['chain_bytes'])
    hm3 = load_hm3(plan['reference'], plan['allele_list'])
    effective_n = 4.0 / (1.0/51643 + 1.0/446273)
    if effective_n != plan['assumed_effective_N']:
        raise RuntimeError('ASSUMED_CONSTANT_N_CHANGED')
    dispositions, diagnostics, seen = Counter(), Counter(), set()
    raw_digest, out_digest, identity_digest = hashlib.sha256(), hashlib.sha256(), hashlib.sha256()
    max_relative = None
    source_rows = 0
    with gzip.open(source, 'rb') as input_handle, output.open('xb') as raw_output:
        with gzip.GzipFile(fileobj=raw_output, mode='wb', filename='', mtime=0, compresslevel=1) as dest:
            first = input_handle.readline()
            raw_digest.update(first)
            if first.decode('utf-8').rstrip('\r\n').split('\t') != HEADER:
                raise RuntimeError('LITERAL_R13_HEADER_CHANGED')
            out_header = b'SNP\tA1\tA2\tZ\tN\n'
            dest.write(out_header); out_digest.update(out_header)
            for line in input_handle:
                raw_digest.update(line)
                source_rows += 1
                values = line.decode('utf-8').rstrip('\r\n').split('\t')
                if len(values) != len(HEADER):
                    raise RuntimeError('R13_SOURCE_ROW_WIDTH_CHANGED')
                disposition, result = qualify_row(dict(zip(HEADER, values)), hm3, chain, seen)
                dispositions[disposition] += 1
                if result is None:
                    continue
                p = result['supplied_p']
                if p is None or not 0 <= p <= 1:
                    diagnostics['supplied_P_missing_or_invalid'] += 1
                elif p == 0:
                    diagnostics['supplied_P_literal_zero'] += 1
                else:
                    computed = math.erfc(abs(result['source_z'])/math.sqrt(2))
                    relative = abs(p-computed)/max(p, computed, 1e-300)
                    diagnostics['supplied_P_finite_positive'] += 1
                    diagnostics['P_Z_relative_difference_gt_0_1'] += int(relative > 0.1)
                    max_relative = relative if max_relative is None else max(max_relative, relative)
                diagnostics['retained_reverse_chain'] += int(result['reverse_chain'])
                seen.add(result['SNP'])
                identity = '\t'.join(result[k] for k in ('SNP','A1','A2')) + '\n'
                identity_digest.update(identity.encode())
                rendered = (identity.rstrip('\n') + '\t' + format(result['Z'], '.12g') + '\t' + format(effective_n, '.12g') + '\n').encode()
                dest.write(rendered); out_digest.update(rendered)
        raw_output.flush()
        import os
        os.fsync(raw_output.fileno())
    if sum(dispositions.values()) != source_rows or len(seen) != dispositions['retained']:
        raise RuntimeError('SOURCE_ROW_ACCOUNTING_FAILED')
    # Verify the complete derivative independently through gzip CRC/EOF.
    derivative_digest, derivative_rows = hashlib.sha256(), 0
    with gzip.open(output, 'rb') as f:
        for line in f:
            derivative_digest.update(line)
            derivative_rows += 1
    if derivative_digest.hexdigest() != out_digest.hexdigest() or derivative_rows != dispositions['retained'] + 1:
        raise RuntimeError('DERIVATIVE_FULL_STREAM_VERIFICATION_FAILED')
    after = body_hashes(source)
    check_hashes(plan['dependencies_sha256'])
    if after != before or sha(plan_path) != plan_hash:
        raise RuntimeError('SOURCE_OR_PLAN_CHANGED_DURING_PREPROCESSING')
    result = dict(schema='qualified_finngen_insomnia_preprocessing_v1', completed_utc=utc(), plan_sha256=plan_hash,
                  scope=plan['scope'], source=str(source), source_before=before, source_after=after,
                  full_source_gzip_CRC_and_EOF_verified=True, source_decompressed_sha256=raw_digest.hexdigest(),
                  literal_header=HEADER, source_rows=source_rows, terminal_dispositions=dict(dispositions),
                  diagnostic_counts=dict(diagnostics), maximum_finite_positive_P_Z_relative_difference=max_relative,
                  P_Z_diagnostics_used_for_filtering=False, reference_nonambiguous_snps=len(hm3),
                  chain_metadata=chain_metadata, derivative=str(output), derivative_sha256=sha(output),
                  derivative_decompressed_sha256=out_digest.hexdigest(), ordered_SNP_allele_sha256=identity_digest.hexdigest(),
                  derivative_full_gzip_CRC_and_EOF_verified=True, derivative_rows=dispositions['retained'],
                  assumed_effective_N=effective_n, verified_per_variant_N=False, verified_per_variant_INFO=False,
                  all_retained_rows_unique_coordinate_and_allele_compatible=True, source_order_preserved=True,
                  retained_minimum_diagnostic_pass=dispositions['retained'] >= 700_000,
                  scientific_source_admitted=False, independent_replication_established=False,
                  status='QUALIFIED_PREPROCESSING_DIAGNOSTIC_PASS' if dispositions['retained'] >= 700_000 else 'FAILED_RETAINED_VARIANT_MINIMUM_PRESERVED')
    write_new(receipt, result)
    print(json.dumps(dict(status=result['status'], source_rows=source_rows, retained_rows=result['derivative_rows'])), flush=True)
    if not result['retained_minimum_diagnostic_pass']:
        raise RuntimeError('HISTORICAL_MINIMUM_RETAINED_VARIANTS_FAILED_NO_FIT')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--plan-sha256', required=True)
    parser.add_argument('--inherited-heavy-lock-fd', type=int, required=True)
    args = parser.parse_args()
    preprocess(args.plan, args.plan_sha256, args.inherited_heavy_lock_fd)


if __name__ == '__main__':
    main()

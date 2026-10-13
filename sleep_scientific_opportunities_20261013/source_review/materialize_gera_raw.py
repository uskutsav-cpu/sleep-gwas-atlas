#!/usr/bin/env python3
"""Exact GERA additive raw37 -> reference-oriented five-column input; no fit.

Implementation only, pending independent admission. SQLite holds reference and
candidate metadata on SSD, never full raw rows. Source P is diagnostic only.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from fractions import Fraction
import gzip
import hashlib
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import sqlite3
import sys

ROOT = Path(__file__).resolve().parents[2]
SSD = Path('/Volumes/Extreme SSD')
HELPER = ROOT / 'discovery_extension/scripts/47_stream_replication_sources.py'
HELPER_INTERFACE_SHA = 'aee2ef31ca1d6fb7c0e956c14e7553db39993f9878981e2a0117658c78e8d256'
REFERENCE = SSD / 'sleep-unified-research-v1/recovery-2026-10-08/native_inputs/discovery_extension/data/reference/panukbb_hm3_variant_reference.tsv.gz'
ALLELES = SSD / 'Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas/work/track_b_completion/local_dependency_copies/ref/eur_w_ld_chr/w_hm3.snplist'
REFERENCE_SHA = 'e6e4814d99a1eff91875014fe70c61f760182efdc1711ba6e155963cf5aa11f8'
ALLELES_SHA = 'ec73fca0b696e8beba465b51e52911676fcade375bcc9475d99c0ec30509d3ed'
SOURCE_URL = 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90086001-GCST90087000/GCST90086097/GCST90086097_buildGRCh37.tsv'
SOURCE_BYTES = 2019290168
SOURCE_MD5 = '6dfa12cf01f443e1115d25b83a77d703'
SCIENTIFIC_PROTOCOL_SHA = 'df1372bb01c683dcaae2a74e6b59ef44b69d4ed9b73ddc32a4908248deb962ed'
HEADER = (
    'chromosome', 'base_pair_location', 'rs_id_all', 'info_all',
    'other_allele', 'effect_allele', 'all_maf', 'all_total',
    'standard_error', 'beta', 'p_value', 'freq_effect_allele',
    'chr_X_freq_alleleB_all', 'chr_X_freq_alleleB_males',
    'chr_X_freq_alleleB_females', 'chr_X_all_total',
    'chr_X_frequentist_add_pvalue_males',
    'chr_X_frequentist_add_beta_1.genotype.sex.1_males',
    'chr_X_frequentist_add_se_1.genotype.sex.1_males',
    'chr_X_all_total_males', 'chr_X_frequentist_add_pvalue_females',
    'chr_X_frequentist_add_beta_1.genotype.sex.2_females',
    'chr_X_frequentist_add_se_1.genotype.sex.2_females', 'chr_X_all_total_females',
)
VALID = {'A', 'C', 'G', 'T'}
PALINDROMES = {frozenset(('A', 'T')), frozenset(('C', 'G'))}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()


def import_orientation():
    """Import the existing function without invoking its worker or creating pyc."""
    text = HELPER.read_text()
    interface = ''.join(text.splitlines(keepends=True)[57:93]).encode()
    if hashlib.sha256(interface).hexdigest() != HELPER_INTERFACE_SHA:
        raise ValueError('ORIENTATION_INTERFACE_CHANGED_REVIEW_REQUIRED')
    old_bytecode = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(HELPER.parent))
    try:
        spec = importlib.util.spec_from_file_location('gera_reviewed_orientation', HELPER)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        if (module.VALID != VALID or module.AMBIGUOUS != PALINDROMES or
                module.COMPLEMENT != str.maketrans('ACGT','TGCA')):
            raise ValueError('ORIENTATION_GLOBAL_CONSTANTS_CHANGED_REVIEW_REQUIRED')
        return module.orientation
    finally:
        sys.path.pop(0)
        sys.dont_write_bytecode = old_bytecode


orientation = import_orientation()


def tsv_rows(path, required, compressed=False):
    opener = gzip.open if compressed else Path.open
    with opener(path, 'rt', encoding='utf-8', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames or []
        if len(set(fields)) != len(fields) or not set(required).issubset(fields):
            raise ValueError('REFERENCE_HEADER_INVALID')
        for row in reader:
            if None in row or any(row[k] is None for k in fields):
                raise ValueError('REFERENCE_ROW_MALFORMED')
            yield row


def load_reference(db, reference, alleles, minimum=900000):
    """Strict old load_hm3 concepts; every compatible coordinate candidate kept."""
    db.executescript('''
        CREATE TABLE allele_definitions(snp TEXT, signature TEXT, a1 TEXT, a2 TEXT,
            PRIMARY KEY(snp,signature));
        CREATE TABLE reference_definitions(snp TEXT, signature TEXT, chr INTEGER, bp INTEGER,
            PRIMARY KEY(snp,signature));
        CREATE TABLE refs(snp TEXT PRIMARY KEY, a1 TEXT, a2 TEXT, chr INTEGER, bp INTEGER);
    ''')
    counts = Counter()
    for row in tsv_rows(alleles, ('SNP', 'A1', 'A2')):
        if not re.fullmatch(r'rs[1-9][0-9]*', row['SNP']):
            raise ValueError('REFERENCE_NOT_AUTHENTICATED_RSID')
        # Exact parsed-field equality, including any additional fields. A JSON
        # signature is the full definition, not a collision-prone hash proof.
        signature = json.dumps(row, sort_keys=True, separators=(',',':'))
        inserted = db.execute('INSERT OR IGNORE INTO allele_definitions VALUES(?,?,?,?)',
                             (row['SNP'], signature, row['A1'].upper(), row['A2'].upper())).rowcount
        if not inserted:
            counts['identical_allele_definition_copies'] += 1
        counts['allele_list_rows'] += 1
    db.commit()
    for row in tsv_rows(reference, ('SNP', 'CHR', 'BP'), compressed=True):
        counts['reference_rows'] += 1
        if db.execute('SELECT 1 FROM allele_definitions WHERE snp=? LIMIT 1', (row['SNP'],)).fetchone() is None:
            counts['reference_not_in_allele_list'] += 1
            continue
        chromosome, position = int(row['CHR']), int(row['BP'])
        signature = json.dumps(row, sort_keys=True, separators=(',',':'))
        inserted = db.execute('INSERT OR IGNORE INTO reference_definitions VALUES(?,?,?,?)',
                             (row['SNP'], signature, chromosome, position)).rowcount
        if not inserted:
            counts['identical_reference_definition_copies'] += 1
        if counts['reference_rows'] % 10000 == 0:
            db.commit()
    db.executescript('''
        CREATE TABLE tainted_snps(snp TEXT PRIMARY KEY);
        INSERT INTO tainted_snps SELECT snp FROM allele_definitions GROUP BY snp HAVING count(*)>1;
        INSERT OR IGNORE INTO tainted_snps
            SELECT snp FROM reference_definitions GROUP BY snp HAVING count(*)>1;
        CREATE TABLE tainted_coordinates AS SELECT DISTINCT r.chr,r.bp
            FROM reference_definitions r JOIN tainted_snps t ON r.snp=t.snp;
        CREATE UNIQUE INDEX tainted_coordinate ON tainted_coordinates(chr,bp);
    ''')
    counts['conflicting_allele_SNP_groups'] = db.execute('''SELECT count(*) FROM
        (SELECT snp FROM allele_definitions GROUP BY snp HAVING count(*)>1)''').fetchone()[0]
    counts['conflicting_reference_SNP_groups'] = db.execute('''SELECT count(*) FROM
        (SELECT snp FROM reference_definitions GROUP BY snp HAVING count(*)>1)''').fetchone()[0]
    counts['tainted_SNPs'] = db.execute('SELECT count(*) FROM tainted_snps').fetchone()[0]
    counts['tainted_coordinate_groups'] = db.execute('SELECT count(*) FROM tainted_coordinates').fetchone()[0]
    for snp, a1, a2, chromosome, position in db.execute('''SELECT r.snp,a.a1,a.a2,r.chr,r.bp
            FROM reference_definitions r JOIN allele_definitions a ON r.snp=a.snp
            WHERE NOT EXISTS (SELECT 1 FROM tainted_snps t WHERE t.snp=r.snp)'''):
        if (chromosome not in range(1, 23) or position <= 0 or a1 not in VALID
                or a2 not in VALID or a1 == a2 or frozenset((a1,a2)) in PALINDROMES):
            counts['reference_ineligible_coordinate_or_alleles'] += 1
            continue
        db.execute('INSERT INTO refs VALUES(?,?,?,?,?)',
                   (snp, a1, a2, chromosome, position))
        counts['eligible_reference_rsids'] += 1
    db.commit()
    db.execute('CREATE INDEX reference_coordinate ON refs(chr,bp)')
    counts['coordinate_collision_groups'] = db.execute('''SELECT count(*) FROM
        (SELECT chr,bp FROM refs GROUP BY chr,bp HAVING count(*)>1)''').fetchone()[0]
    if counts['eligible_reference_rsids'] < minimum:
        raise ValueError('UNEXPECTEDLY_SMALL_NONAMBIGUOUS_REFERENCE')
    return dict(counts)


def finite(token):
    try:
        value = float(token)
        return value if math.isfinite(value) else None
    except (ValueError, TypeError):
        return None


def total_n(token):
    """Source decimal integer test; float rounding must never admit a fraction."""
    if token.strip() in ('', '.', 'NA', 'N/A', 'null'):
        return None, 'N_missing'
    try:
        value = Decimal(token)
    except InvalidOperation:
        return None, 'N_invalid_decimal'
    if not value.is_finite():
        return None, 'N_nonfinite'
    if value != value.to_integral_value():
        return None, 'N_fractional'
    if value <= 0:
        return None, 'N_nonpositive'
    if value > 56637:
        return None, 'N_above_56637'
    return int(value), None


def classify(row, diagnostics, n_hist):
    """Return structural disposition separately from deferred statistical filter.

    Every uniquely mappable structural row enters duplicate grouping, including
    rows whose beta/SE, INFO/MAF or N will fail after grouping.
    """
    info, maf = (finite(row[k]) for k in ('info_all', 'all_maf'))
    n, n_reason = total_n(row['all_total'])
    beta, se, p = (finite(row[k]) for k in ('beta', 'standard_error', 'p_value'))
    for name, value in [('INFO', info), ('MAF', maf)]:
        if value is None:
            diagnostics['missing_or_nonfinite_'+name] += 1
        else:
            low, high = 'min_finite_'+name, 'max_finite_'+name
            diagnostics[low] = value if low not in diagnostics else min(value, diagnostics[low])
            diagnostics[high] = value if high not in diagnostics else max(value, diagnostics[high])
    if info is not None and not 0 <= info <= 1:
        diagnostics['INFO_outside_0_1'] += 1
    if maf is not None and not 0 <= maf <= 0.5:
        diagnostics['MAF_outside_0_0.5'] += 1
    if info is not None and info <= 0.9:
        diagnostics['INFO_not_gt_0.9'] += 1
    if maf is not None and not 0.01 < maf <= 0.5:
        diagnostics['MAF_not_in_0.01_0.5'] += 1
    if beta is None:
        diagnostics['beta_missing_or_nonfinite'] += 1
    if se is None:
        diagnostics['SE_missing_or_nonfinite'] += 1
    elif se <= 0:
        diagnostics['SE_not_positive'] += 1
    if beta is not None and se is not None and se > 0 and not math.isfinite(beta/se):
        diagnostics['nonfinite_supplied_Z'] += 1
    if n_reason:
        diagnostics[n_reason] += 1
    else:
        n_hist[n] += 1
    if row['rs_id_all'] in ('', '.', 'NA'):
        diagnostics['source_ID_missing'] += 1
    if p is None or not 0 <= p <= 1:
        diagnostics['source_P_missing_or_invalid'] += 1
    elif beta is not None and se is not None and se > 0:
        z = beta / se
        if math.isfinite(z):
            wald = math.erfc(abs(z) / math.sqrt(2))
            diagnostics['source_P_Wald_diagnostic_rows'] += 1
            diagnostics['max_abs_source_P_minus_Wald_P'] = max(
                diagnostics['max_abs_source_P_minus_Wald_P'], abs(p-wald))
    try:
        chromosome = int(row['chromosome'])
        position = int(row['base_pair_location'])
    except ValueError:
        return 'invalid_integer_coordinate', None
    if chromosome not in range(1, 23):
        return 'not_autosome', None
    if position <= 0:
        return 'invalid_positive_1based_position', None
    if chromosome == 6 and 25000000 <= position <= 34000000:
        return 'extended_MHC37', None
    effect, other = row['effect_allele'].upper(), row['other_allele'].upper()
    if effect not in VALID or other not in VALID or effect == other:
        return 'not_distinct_single_base_ACGT', None
    if frozenset((effect, other)) in PALINDROMES:
        return 'palindrome', None
    z = beta/se if beta is not None and se is not None and se > 0 else None
    if beta is None:
        reason = 'beta_missing_or_nonfinite'
    elif se is None:
        reason = 'SE_missing_or_nonfinite'
    elif se <= 0:
        reason = 'SE_not_positive'
    elif not math.isfinite(z):
        reason = 'nonfinite_supplied_Z'
    elif info is None:
        reason = 'INFO_missing_or_nonfinite'
    elif info <= 0.9:
        reason = 'INFO_not_gt_0.9'
    elif maf is None:
        reason = 'MAF_missing_or_nonfinite'
    elif not 0.01 < maf <= 0.5:
        reason = 'MAF_not_in_0.01_0.5'
    else:
        reason = n_reason
    return None, (chromosome, position, effect, other, z, n, reason)


def scan_source(db, source, identity):
    db.executescript('''
        CREATE TABLE candidates(snp TEXT, a1 TEXT, a2 TEXT, z REAL, n INTEGER,
            seq INTEGER, offset INTEGER, length INTEGER, digest TEXT, copies INTEGER,
            exclusion TEXT);
        CREATE INDEX candidate_snp_digest ON candidates(snp,digest);
    ''')
    excluded, diagnostics, n_hist, orientations = Counter(), Counter(), Counter(), Counter()
    input_rows, duplicates = 0, 0
    h_md5, h_sha = hashlib.md5(), hashlib.sha256()
    with source.open('rb') as f, source.open('rb') as duplicate_reader:
        header = f.readline()
        if header.rstrip(b'\r\n').decode('utf-8') != '\t'.join(HEADER):
            raise ValueError('LITERAL_GERA_HEADER_MISMATCH')
        h_md5.update(header)
        h_sha.update(header)
        while True:
            offset = f.tell()
            raw = f.readline()
            if not raw:
                break
            h_md5.update(raw)
            h_sha.update(raw)
            input_rows += 1
            try:
                fields = next(csv.reader([raw.decode('utf-8')], delimiter='\t', strict=True))
            except (UnicodeDecodeError, csv.Error, StopIteration):
                excluded['malformed_record'] += 1
                continue
            if len(fields) != len(HEADER):
                excluded['malformed_record'] += 1
                continue
            row = dict(zip(HEADER, fields))
            reason, values = classify(row, diagnostics, n_hist)
            if reason:
                excluded[reason] += 1
                continue
            chromosome, position, effect, other, z, n, deferred_exclusion = values
            if db.execute('SELECT 1 FROM tainted_coordinates WHERE chr=? AND bp=?',
                          (chromosome, position)).fetchone() is not None:
                excluded['tainted_reference_coordinate'] += 1
                continue
            compatible = []
            coordinate_present = False
            for snp, a1, a2 in db.execute('SELECT snp,a1,a2 FROM refs WHERE chr=? AND bp=?', (chromosome, position)):
                coordinate_present = True
                sign = orientation(effect, other, a1, a2)
                if sign is not None:
                    compatible.append((snp, a1, a2, sign))
            if len(compatible) != 1:
                category = ('multiple_compatible_rsids' if compatible else
                            'reference_allele_mismatch' if coordinate_present else 'outside_reference_position')
                excluded[category] += 1
                continue
            snp, a1, a2, sign = compatible[0]
            pair = (effect, other)
            category = ('same' if pair == (a1, a2) else 'swap' if pair == (a2, a1)
                        else 'complement' if sign == 1 else 'complement_swap')
            orientations[category] += 1
            fingerprint = hashlib.sha256(raw).hexdigest()
            identical = None
            for rid, old_offset, old_length in db.execute('''SELECT rowid,offset,length FROM
                    candidates WHERE snp=? AND digest=?''', (snp, fingerprint)):
                duplicate_reader.seek(old_offset)
                if old_length == len(raw) and duplicate_reader.read(old_length) == raw:
                    identical = rid
                    break
            if identical is not None:
                db.execute('UPDATE candidates SET copies=copies+1 WHERE rowid=?', (identical,))
                duplicates += 1
            else:
                db.execute('INSERT INTO candidates VALUES(?,?,?,?,?,?,?,?,?,?,?)',
                           (snp, a1, a2, sign*z if z is not None else None, n, input_rows,
                            offset, len(raw), fingerprint, 1, deferred_exclusion))
            if input_rows % 10000 == 0:
                db.commit()
        actual_bytes = f.tell()
    db.commit()
    actual = {'bytes': actual_bytes, 'md5': h_md5.hexdigest(), 'sha256': h_sha.hexdigest()}
    if any(actual[k] != identity[k] for k in actual):
        raise ValueError('PROCESSED_SOURCE_IDENTITY_MISMATCH_NO_OUTPUT_ADMITTED')
    return {'input_rows': input_rows, 'exclusive_exclusions': dict(excluded),
            'duplicate_copies': duplicates, 'diagnostic_flags': dict(diagnostics),
            'orientation_counts_mapped_input_rows': dict(orientations),
            'source_header_sha256': hashlib.sha256(header).hexdigest(),
            'source_identity': actual, 'valid_total_N_histogram': n_hist}


def r7_histogram(histogram, p):
    """Integer-valued N, exact rational R7 interpolation, without expanding rows."""
    count = sum(histogram.values())
    if not count:
        return None
    rank = p*(count-1)
    low, remainder = divmod(rank.numerator, rank.denominator)
    high = low+(remainder != 0)
    cumulative, lo, hi = 0, None, None
    for n, copies in sorted(histogram.items()):
        cumulative += copies
        if lo is None and cumulative > low:
            lo = n
        if cumulative > high:
            hi = n
            break
    return Fraction(lo)+(Fraction(remainder, rank.denominator)*(hi-lo))


def emit_unique(db, output):
    db.execute('''CREATE TABLE groups AS SELECT snp,count(*) AS variants,sum(copies) AS
                  source_rows FROM candidates GROUP BY snp''')
    db.execute('CREATE UNIQUE INDEX group_snp ON groups(snp)')
    conflicts, conflict_records, conflict_source_rows = db.execute('''SELECT count(*),
        coalesce(sum(variants),0),coalesce(sum(source_rows),0) FROM groups WHERE variants>1''').fetchone()
    post_group_exclusions = dict(db.execute('''SELECT c.exclusion,count(*) FROM candidates c
        JOIN groups g ON c.snp=g.snp WHERE g.variants=1 AND c.exclusion IS NOT NULL GROUP BY c.exclusion'''))
    histogram = Counter(dict(db.execute('''SELECT c.n,count(*) FROM candidates c JOIN groups g
        ON c.snp=g.snp WHERE g.variants=1 AND c.exclusion IS NULL GROUP BY c.n''')))
    population = sum(histogram.values())
    q90 = r7_histogram(histogram, Fraction(9, 10))
    cutoff = None if q90 is None else q90/Fraction(3, 2)
    Q = None if q90 is None else int(10*q90)
    j, remainder = (None, None) if not population else divmod(9*(population-1), 10)
    minimum_n = None if Q is None else (Q+14)//15
    retained, sum_z2, nonfinite_chisq = 0, 0.0, 0
    if population:
        with output.open('xb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0) as gz, io.TextIOWrapper(gz, encoding='utf-8', newline='') as f:
            writer = csv.writer(f, delimiter='\t', lineterminator='\n')
            writer.writerow(('SNP', 'A1', 'A2', 'Z', 'N'))
            if minimum_n is not None:
                for snp, a1, a2, z, n in db.execute('''SELECT c.snp,c.a1,c.a2,c.z,c.n
                        FROM candidates c JOIN groups g ON c.snp=g.snp
                        WHERE g.variants=1 AND c.exclusion IS NULL AND 15*c.n>=? ORDER BY c.seq''', (Q,)):
                    writer.writerow((snp, a1, a2, format(z, '.17g'), n))
                    retained += 1
                    z2 = z*z
                    if math.isfinite(z2) and math.isfinite(sum_z2+z2):
                        sum_z2 += z2
                    else:
                        nonfinite_chisq += 1
    return {'conflicting_groups': conflicts, 'conflicting_unique_records': conflict_records,
            'conflicting_source_rows_diagnostic': conflict_source_rows,
            'post_group_exclusive_exclusions': post_group_exclusions,
            'N_quantile_population_unique_nonconflicting_eligible_mapped': population,
            'R7_exact_m_j_r_Q': {'m': population, 'j': j, 'r': remainder, 'Q': Q},
            'q90_R7': None if q90 is None else float(q90),
            'q90_exact_fraction': None if q90 is None else str(q90),
            'N_cutoff_q90_div_1.5': None if cutoff is None else float(cutoff),
            'minimum_retained_integer_N': minimum_n, 'N_filtered_unique_records': population-retained,
            'N_below_threshold': sum(c for n,c in histogram.items() if 15*n<Q) if Q is not None else 0,
            'N_equal_threshold': sum(c for n,c in histogram.items() if 15*n==Q) if Q is not None else 0,
            'N_above_threshold': sum(c for n,c in histogram.items() if 15*n>Q) if Q is not None else 0,
            'retained_N_min': min((n for n in histogram if 15*n>=Q), default=None) if Q is not None else None,
            'retained_N_max': max((n for n in histogram if 15*n>=Q), default=None) if Q is not None else None,
            'retained_rows': retained, 'nonfinite_chi_square_diagnostic_rows': nonfinite_chisq,
            'mean_chi_square_diagnostic': sum_z2/retained if retained and not nonfinite_chisq else None}


def materialize(source, reference, alleles, directory, identity, minimum_reference=900000):
    """Called on synthetic fixtures in tests; production CLI is separately gated."""
    with source.open('rb') as f:
        if f.readline().rstrip(b'\r\n').decode('utf-8') != '\t'.join(HEADER):
            raise ValueError('LITERAL_GERA_HEADER_MISMATCH')
    directory.mkdir(mode=0o700, exist_ok=False)
    os.environ['SQLITE_TMPDIR'] = str(directory)
    # This SSD returns SQLITE_READONLY_DBMOVED with the default inode-checking
    # VFS. unix-dotfile retains locking; unix-none (no locks) is not used.
    db = sqlite3.connect((directory/'gera_metadata.sqlite').resolve().as_uri()+'?vfs=unix-dotfile', uri=True)
    try:
        db.execute('PRAGMA temp_store=FILE')
        db.execute('PRAGMA cache_size=-8192')
        ref_counts = load_reference(db, reference, alleles, minimum_reference)
        qc = scan_source(db, source, identity)
        histogram = qc.pop('valid_total_N_histogram')
        qc['valid_source_total_N_summary'] = {
            'rows': sum(histogram.values()), 'distinct_values': len(histogram),
            'min': min(histogram) if histogram else None,
            'max': max(histogram) if histogram else None,
            'R7_quantiles': {str(p): None if not histogram else float(r7_histogram(histogram, p))
                             for p in (Fraction(1,10), Fraction(1,2), Fraction(9,10))}}
        output = directory/'GERA_GCST90086097.supplied_Z.total_N.sumstats.gz'
        qc.update(emit_unique(db, output))
        post_group_exclusions = qc.pop('post_group_exclusive_exclusions')
        exclusions = Counter(qc['exclusive_exclusions'])
        exclusions.update(post_group_exclusions)
        qc['exclusive_exclusions'] = dict(exclusions)
        reconciled = (sum(qc['exclusive_exclusions'].values())+qc['duplicate_copies']+
                      qc['conflicting_unique_records']+qc['N_filtered_unique_records']+qc['retained_rows'])
        if reconciled != qc['input_rows']:
            raise AssertionError('INPUT_DISPOSITION_RECONCILIATION_FAILED')
        qc.update({'schema': 'GERA_raw37_preprocessing_QC/1',
                   'status': ('PREPROCESSING_COMPLETE_PENDING_INDEPENDENT_ACTUAL_OUTPUT_REVIEW'
                              if qc['retained_rows'] else 'HOLD_NO_NATIVE_INPUT'),
                   'input_reconciled': True, 'reconciled_rows': reconciled,
                   'exclusive_precedence': ['malformed_record', 'coordinate/autosome', 'MHC37',
                       'distinct_ACGT', 'palindrome', 'reference_compatibility',
                       'structurally_mapped_duplicate_copies', 'conflicting_unique_records',
                       'beta/SE/Z', 'INFO', 'MAF', 'N', 'N_missingness', 'retained'],
                   'duplicate_policy': 'Group ALL structurally valid uniquely mapped rows BEFORE statistical filters. SHA matches require source-offset full-byte proof. Different records exclude the whole SNP group.',
                   'reference_counts': ref_counts, 'MHC37_inclusive': [6,25000000,34000000],
                   'N_convention': 'PER_VARIANT_ANALYSED_TOTAL_all_total; never Neff or a constant',
                   'Z_convention': 'reference_orientation*(supplied_additive_beta/supplied_SE)',
                   'source_P_used_for_selection_or_Z': False,
                   'reference_provenance': {'path': str(reference),
                       'prior_authenticated_sha256': REFERENCE_SHA if reference.resolve()==REFERENCE.resolve() else 'SYNTHETIC_FIXTURE',
                       'alleles': str(alleles),
                       'alleles_prior_authenticated_sha256': ALLELES_SHA if alleles.resolve()==ALLELES.resolve() else 'SYNTHETIC_FIXTURE',
                       'rehash_performed': False},
                   'reference_duplicate_policy': 'Exact full parsed-field duplicate definitions collapse and count. Different reference or allele-list definitions taint their SNP and ALL associated coordinates; affected source coordinates are held, with no first/last winner.',
                   'SQLite_staging': {'VFS': 'unix-dotfile', 'temp_store': 'FILE', 'cache_KiB': 8192,
                                      'raw_payloads_stored': False, 'scratch_directory': str(directory)},
                   'source_ID_coordinate_consistency': 'UNKNOWN_FORMAT_NOT_DOCUMENTED; generic provider ID description does not specify a parsing grammar; no guessed repair.',
                   'orientation_interface_sha256': HELPER_INTERFACE_SHA,
                   'sumstats_sha256': sha256(output) if output.exists() else None,
                   'sumstats_path': str(output) if output.exists() else None,
                   'native_fits': 0, 'disease_pair_tests': 0, 'replication_claims': 0,
                   'h2_or_power': 'UNKNOWN_NOT_TESTED', 'recorded_utc': datetime.now(timezone.utc).isoformat()})
        db.commit()
        return qc
    finally:
        db.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--acquisition-receipt', type=Path, required=True)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--admission', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    review = json.loads(args.admission.read_text())
    if (review.get('preprocessing_admitted') is not True or
            review.get('reviewed_script_sha256') != sha256(Path(__file__)) or
            review.get('reviewed_protocol_sha256') != sha256(args.protocol) or
            sha256(args.protocol) != SCIENTIFIC_PROTOCOL_SHA):
        raise ValueError('INDEPENDENT_PREPROCESSING_ADMISSION_REQUIRED')
    if (not args.output_dir.resolve().is_relative_to(SSD.resolve()) or
            not args.source.resolve().is_relative_to(SSD.resolve())):
        raise ValueError('SCRATCH_AND_OUTPUT_MUST_BE_ON_SSD')
    if shutil.disk_usage(ROOT).free < 3*1024**3 or shutil.disk_usage(SSD).free < 5*1024**3:
        raise ValueError('UNCHANGED_PREPROCESSING_RESOURCE_FLOOR')
    receipt = json.loads(args.acquisition_receipt.read_text())
    if (receipt.get('status') != 'EXACT_SOURCE_ACQUISITION_PASS' or
            receipt.get('url') != SOURCE_URL or
            Path(receipt.get('path','')).resolve() != args.source.resolve() or
            receipt.get('expected_bytes') != SOURCE_BYTES or receipt.get('actual_bytes') != SOURCE_BYTES or
            receipt.get('actual_md5') != SOURCE_MD5 or receipt.get('expected_md5') != SOURCE_MD5 or
            args.source.stat().st_size != SOURCE_BYTES or
            not re.fullmatch(r'[0-9a-f]{64}', receipt.get('actual_sha256',''))):
        raise ValueError('COMPLETE_EXACT_ACQUISITION_REQUIRED')
    # The acquisition receipt verifies the complete body before preparation.
    # scan_source rehashes every byte to verify the actual processed stream before emission.
    identity = {'bytes': SOURCE_BYTES, 'md5': SOURCE_MD5, 'sha256': receipt['actual_sha256']}
    # Future execution checks actual reference identity before source census.
    # Authoring/tests do not read these files; fixture callers use synthetic inputs.
    if sha256(REFERENCE) != REFERENCE_SHA or sha256(ALLELES) != ALLELES_SHA:
        raise ValueError('AUTHENTICATED_REFERENCE_IDENTITY_REQUIRED')
    os.umask(0o077)
    qc = materialize(args.source, REFERENCE, ALLELES, args.output_dir, identity)
    qc.update({'script_sha256': sha256(Path(__file__)), 'protocol_sha256': sha256(args.protocol),
               'admission_sha256': sha256(args.admission),
               'acquisition_receipt_sha256': sha256(args.acquisition_receipt),
               'accession': 'GCST90086097', 'source_URL': SOURCE_URL})
    qc['reference_provenance']['rehash_performed'] = True
    with (args.output_dir/'GERA_preprocessing.qc.json').open('x') as f:
        json.dump(qc, f, indent=2, sort_keys=True)
        f.write('\n')
    print(json.dumps({'status': qc['status'], 'input_rows': qc['input_rows'],
                      'retained_rows': qc['retained_rows'], 'output_dir': str(args.output_dir)}))


if __name__ == '__main__':
    main()

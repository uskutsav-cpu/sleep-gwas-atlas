"""Identity-only spatial/physical screen for additional measured neuron activity.

This does not change a fit or authorize the screened constructs for training.
Full source components are built before using measurement availability. Every
historical region and the entire reserved study are protected conservatively.
"""
import ast
import bisect
import collections
import csv
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import json
from pathlib import Path
import struct
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/sleep_frailty_ml_neuronal_activity'
CHANNEL = ROOT / 'outputs/sleep_frailty_ml_channels/source_candidates'
sys.path.append(str(ROOT / 'work/stage5_env/lib/python3.13/site-packages'))
from pyliftover import LiftOver


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def rc(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def rows(path, delimiter=','):
    path = Path(path)
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as handle:
        yield from csv.DictReader(handle, delimiter=delimiter)


def integer(value):
    number = Decimal(value)
    if not number.is_finite() or number != number.to_integral_value():
        raise ValueError(value)
    return int(number)


class Components:
    def __init__(self, n):
        self.parent = list(range(n))

    def find(self, i):
        while i != self.parent[i]:
            self.parent[i] = self.parent[self.parent[i]]
            i = self.parent[i]
        return i

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.parent[max(a, b)] = min(a, b)


def main():
    started = time.monotonic()
    original = json.loads((OUT / 'additional_neuronal_activity_inventory_receipt.json').read_text())
    metadata = OUT / Path(original['sources'][0]['path']).name
    assert sha(metadata) == original['sources'][0]['sha256']
    inventory = OUT / 'original_neuronal_activity_identity_inventory.csv.gz'
    assert sha(inventory) == original['inventory_sha256']
    available = {r['name']: r for r in rows(inventory)}
    data = ROOT / 'outputs/sleep_frailty_ml_next/coverage_training_arrays.npz'
    assert sha(data) == original['current_training_arrays_sha256']
    source = list(rows(metadata, '\t'))
    assert len({r['name'] for r in source}) == len(source) == 80215
    comp = Components(len(source))
    canonical = {'chr' + str(i) for i in range(1, 23)} | {'chrX', 'chrY'}
    forward_path = ROOT / 'work/hg19ToHg38.over.chain.gz'
    backward_path = ROOT / 'work/hg38ToHg19.over.chain.gz'
    forward, backward = LiftOver(str(forward_path)), LiftOver(str(backward_path))
    intervals = collections.defaultdict(list)
    mapped = {}
    mapping = collections.Counter()
    physical = {}
    for i, row in enumerate(source):
        sequence = row['sequence'].upper()
        if not set(sequence) - set('ACGT'):
            key = min(sequence, rc(sequence))
            if key in physical:
                comp.union(i, physical[key])
            else:
                physical[key] = i
        try:
            begin, end = integer(row['start']), integer(row['end'])
            chrom = row['chr']
            assert chrom in canonical and 0 <= begin < end
        except (ValueError, InvalidOperation, AssertionError):
            mapping['missing_or_invalid_declared_interval'] += 1
            continue
        if row['ref'] == 'GRCh37':
            a = forward.convert_coordinate(chrom, begin)
            b = forward.convert_coordinate(chrom, end - 1)
            if len(a) != 1 or len(b) != 1 or (a[0][0], a[0][2]) != (b[0][0], b[0][2]):
                mapping['nonunique_or_inconsistent_GRCh37_conversion'] += 1
                continue
            c, p, strand, _ = a[0]
            _, q, _, _ = b[0]
            ba, bb = backward.convert_coordinate(c, p), backward.convert_coordinate(c, q)
            if (len(ba) != 1 or len(bb) != 1 or ba[0][:3] != (chrom, begin, strand)
                    or bb[0][:3] != (chrom, end - 1, strand) or abs(q - p) + 1 != end - begin):
                mapping['nonreciprocal_or_length_changed_GRCh37_conversion'] += 1
                continue
            chrom, begin, end = c, min(int(p), int(q)), max(int(p), int(q)) + 1
            mapping['reciprocal_length_preserving_GRCh37_conversion'] += 1
        elif row['ref'] != 'GRCh38':
            mapping['unsupported_declared_build'] += 1
            continue
        else:
            mapping['declared_GRCh38_interval'] += 1
        if chrom not in canonical:
            mapping['noncanonical_mapped_chromosome'] += 1
            continue
        mapped[i] = (chrom, begin, end)
        intervals[chrom].append((begin, end, i))
    # All categories, missing measurements and reference builds participate.
    for chrom, values in intervals.items():
        previous, farthest = None, -1
        for begin, end, i in sorted(values):
            if previous is not None and begin <= farthest + 2000:
                comp.union(previous, i)
            else:
                previous = i
            farthest = max(farthest, end)
    print('COMPONENTS BUILT', len(source), len({comp.find(i) for i in range(len(source))}), flush=True)
    audit = json.loads((CHANNEL / 'Ulirsch_physical_and_source_count_audit.json').read_text())
    history = dict(audit['historical_files_sha256'])
    extra = ROOT / 'outputs/sleep_frailty_ml_next/OL27_development_sequence_bank.csv.gz'
    receipt = next(x for x in original['physical_history_sources'] if x['path'] == str(extra))
    assert sha(extra) == receipt['sha256']
    history[str(extra)] = receipt['sha256']
    past = collections.defaultdict(set)
    history_profile = []
    for filename, expected in history.items():
        assert sha(filename) == expected
        count = 0
        for row in rows(filename):
            past[row['chromosome']].add(int(row['position_hg38']) - 1)
            count += 1
        history_profile.append({'path': filename, 'sha256': expected, 'rows': count})
    past = {c: sorted(values) for c, values in past.items()}
    reserved_path = CHANNEL / 'GSE87711_design.txt.gz'
    retrieval = json.loads((CHANNEL / 'Ulirsch_original_retrieval_receipt.json').read_text())
    expected = next(x['sha256'] for x in retrieval['files'] if Path(x['path']).name == reserved_path.name)
    assert sha(reserved_path) == expected
    reserved = collections.defaultdict(set)
    unmapped = 0
    with gzip.open(reserved_path, 'rt') as handle:
        variants = {identifier.split('_')[0] for identifier, sequence in csv.reader(handle, delimiter='\t')}
    for variant in variants:
        chrom, p, ref, alt = variant.split(':')
        positions = forward.convert_coordinate('chr' + chrom.removeprefix('chr'), int(p) - 1)
        if not positions:
            unmapped += 1
        # Conservatively protect every mapping, without selecting by outcomes.
        for c, q, strand, score in positions:
            reserved[c].add(int(q))
    reserved = {c: sorted(values) for c, values in reserved.items()}
    bad_history, bad_reserved = set(), set()
    for i, (chrom, begin, end) in mapped.items():
        for points, bad, envelope in [(past, bad_history, 2001), (reserved, bad_reserved, 2146)]:
            values = points.get(chrom, [])
            j = bisect.bisect_left(values, begin - envelope)
            if j < len(values) and values[j] <= end + envelope:
                bad.add(comp.find(i))
    # Use the existing independently verified official-genome reader, without
    # importing its model/environment setup or executing its helper functions.
    reader_path = ROOT / 'work/sleep_frailty_ml_v7/reference.py'
    tree = ast.parse(reader_path.read_text())
    cls = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == 'TwoBit')
    namespace = {'struct': struct}
    exec(compile(ast.Module(body=[cls], type_ignores=[]), str(reader_path), 'exec'), namespace)
    genome_path = Path('/Volumes/Extreme SSD/codex-sleep-frailty-work/sleep_frailty_ml_v6/sources/official_hg38.2bit')
    assert sha(genome_path) == audit['genome_sha256']
    genome = namespace['TwoBit'](genome_path)
    records = []
    counts = collections.Counter()
    for i, row in enumerate(source):
        info = available.get(row['name'])
        if info is None or info['history_class'] != 'new_physical_sequence':
            continue
        root = comp.find(i)
        reasons = []
        if root in bad_history:
            reasons.append('whole_source_component_near_any_historical_region')
        if root in bad_reserved:
            reasons.append('whole_source_component_near_entire_reserved_study')
        if i not in mapped:
            reasons.append('no_usable_declared_or_mapped_interval')
        if row['category'] != 'element':
            reasons.append('outside_reference_genomic_element_auxiliary_candidate_scope')
        # This first candidate tranche avoids inferred alternate edits and
        # unresolved build labels. These rows are retained in the audit table.
        if row['ref'] != 'GRCh38':
            reasons.append('not_originally_declared_GRCh38')
        matched = False
        if i in mapped and row['category'] == 'element' and row['ref'] == 'GRCh38':
            chrom, begin, end = mapped[i]
            seq = row['sequence'].upper()
            if end - begin == len(seq):
                truth = genome.sequence(chrom, begin, end)
                matched = seq == truth or rc(seq) == truth
            if not matched:
                reasons.append('full_actual_construct_does_not_match_declared_GRCh38_interval')
        for reason in reasons:
            counts[reason] += 1
        if not reasons:
            counts['passes_this_spatial_and_reference_screen'] += 1
        records.append({**info, 'source_component': root,
            'mapped_GRCh38_interval': ':'.join(map(str, mapped[i])) if i in mapped else '',
            'full_actual_reference_match': matched, 'passes_spatial_and_reference_screen': not reasons,
            'screen_exclusions': '|'.join(reasons), 'training_eligible': False,
            'remaining': 'Shorter shared-DNA core/component audit and duplicate-activity handling; independent guard replay; future experiment plan.'})
    target = OUT / 'additional_activity_spatial_and_reference_screen.csv.gz'
    with gzip.open(target, 'wt') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    result = {'passed_identity_and_reference_screen': True, 'source_rows_before_measurement_filter': len(source),
        'all_source_spatial_and_full_DNA_components': len({comp.find(i) for i in range(len(source))}),
        'component_policy': 'Transitive intervals within 2kb plus exact whole construct DNA in either orientation, across the entire source before coefficient availability. All historical regions, including training, are protected conservatively. Shorter shared-core joins remain pending.',
        'metadata_sha256': sha(metadata), 'source_inventory_sha256': sha(inventory), 'mapping': dict(mapping),
        'historical_sources': history_profile, 'reserved_entire_source_variants': len(variants),
        'reserved_variants_without_GRCh38_mapping': unmapped, 'reserved_original_design_sha256': expected,
        'reserved_genomic_guard': 'Whole-source exclusion if any member is within 2146bp of any mapped original reserved-study variant; includes 145bp construct envelope and a one-base boundary allowance.',
        'new_literal_construct_rows': len(records), 'overlapping_exclusion_counts': dict(counts),
        'potential_reference_element_distinct_physical_sequences': len({r['physical_sequence_sha256_both_orientations'] for r in records if r['passes_spatial_and_reference_screen']}),
        'potential_reference_element_components': len({r['source_component'] for r in records if r['passes_spatial_and_reference_screen']}),
        'coefficient_size_or_P_used_for_screening': False, 'reserved_measurement_files_read': False,
        'genome_sha256': audit['genome_sha256'], 'reference_reader_sha256': sha(reader_path),
        'chain_sha256': {p.name: sha(p) for p in [forward_path, backward_path]},
        'screen_table_sha256': sha(target), 'training_authorized_by_this_screen': False,
        'independent_validation_or_new_model_success_established': False, 'current_comparison_changed': False,
        'remaining': ['Shorter shared-DNA core/component joins against all source constructs, every historical DNA source and the whole reserved study',
                      'Duplicated constructs with conflicting fitted activity coefficients', 'Independent guard replay',
                      'Prespecified additional-activity experiment after current finite comparison, if scientifically justified'],
        'code_sha256': sha(Path(__file__)), 'elapsed_seconds': time.monotonic() - started}
    assert sha(data) == original['current_training_arrays_sha256']
    (OUT / 'additional_activity_spatial_and_reference_receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()

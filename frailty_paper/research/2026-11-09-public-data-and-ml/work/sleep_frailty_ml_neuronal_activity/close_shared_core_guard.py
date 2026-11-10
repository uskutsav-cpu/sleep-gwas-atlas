"""Close exact shared-145bp components and protect all prior/reserved DNA.

No phenotype, reserved effect estimate, p-value or model output is read.
This is preparation for a separately specified future activity experiment.
"""
import bisect
import collections
import csv
from decimal import Decimal
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/sleep_frailty_ml_neuronal_activity'
CHANNEL = ROOT / 'outputs/sleep_frailty_ml_channels/source_candidates'
sys.path.append(str(ROOT / 'work/stage5_env/lib/python3.13/site-packages'))
from pyliftover import LiftOver


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def rc(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def rows(path, delimiter=','):
    path = Path(path)
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        yield from csv.DictReader(f, delimiter=delimiter)


def integral(s):
    d = Decimal(s)
    if not d.is_finite() or d != d.to_integral_value():
        raise ValueError(s)
    return int(d)


def independent_initial_components(source, forward, backward):
    """Interval blocks + graph traversal, rather than the original row unions."""
    canonical = {'chr' + str(i) for i in range(1, 23)} | {'chrX', 'chrY'}
    mapped, by_chrom = {}, collections.defaultdict(list)
    for i, r in enumerate(source):
        try:
            a, b, c = integral(r['start']), integral(r['end']), r['chr']
        except (ValueError, ArithmeticError):
            continue
        if c not in canonical or not 0 <= a < b:
            continue
        if r['ref'] == 'GRCh37':
            left, right = forward.convert_coordinate(c, a), forward.convert_coordinate(c, b - 1)
            if len(left) != 1 or len(right) != 1:
                continue
            l, rr = left[0], right[0]
            if l[0] != rr[0] or l[2] != rr[2] or abs(rr[1] - l[1]) != b - a - 1:
                continue
            ba, bb = backward.convert_coordinate(l[0], l[1]), backward.convert_coordinate(rr[0], rr[1])
            if len(ba) != 1 or len(bb) != 1 or ba[0][:3] != (c, a, l[2]) or bb[0][:3] != (c, b - 1, rr[2]):
                continue
            c, a, b = l[0], int(min(l[1], rr[1])), int(max(l[1], rr[1])) + 1
        elif r['ref'] != 'GRCh38':
            continue
        if c in canonical:
            mapped[i] = (c, a, b)
            by_chrom[c].append((a, b, i))
    blocks, assignment = [], {}
    for values in by_chrom.values():
        current, maximum = None, -1
        for a, b, i in sorted(values):
            if current is None or a > maximum + 2000:
                blocks.append([])
                current = len(blocks) - 1
            blocks[current].append(i)
            assignment[i] = current
            maximum = max(maximum, b)
    for i in range(len(source)):
        if i not in assignment:
            assignment[i] = len(blocks)
            blocks.append([i])
    graph = [set() for block in blocks]
    seen = {}
    for i, r in enumerate(source):
        s = r['sequence'].upper()
        if set(s) - set('ACGT'):
            continue
        key = min(s, rc(s))
        block = assignment[i]
        other = seen.setdefault(key, block)
        if block != other:
            graph[block].add(other)
            graph[other].add(block)
    labels = [None] * len(source)
    visited = set()
    for initial in range(len(blocks)):
        if initial in visited:
            continue
        stack, members = [initial], []
        visited.add(initial)
        while stack:
            block = stack.pop()
            members.extend(blocks[block])
            for neighbor in graph[block] - visited:
                visited.add(neighbor)
                stack.append(neighbor)
        root = min(members)
        for i in members:
            labels[i] = root
    assert all(value is not None for value in labels)
    return labels, mapped


class Join:
    def __init__(self, labels):
        self.parent = labels[:]
        self.merges = 0

    def find(self, i):
        while self.parent[i] != i:
            self.parent[i] = self.parent[self.parent[i]]
            i = self.parent[i]
        return i

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.parent[max(a, b)] = min(a, b)
            self.merges += 1


def cores(sequence):
    reverse = rc(sequence)
    n = len(sequence)
    for j in range(n - 144):
        yield min(sequence[j:j + 145], reverse[n - j - 145:n - j])


def build_lookup(unique_source, join, candidate_ids, protected=()):
    wanted = {join.find(i) for i in candidate_ids} - set(protected)
    lookup = {}
    for sequence, i in unique_source.items():
        if join.find(i) not in wanted:
            continue
        for key in cores(sequence):
            other = lookup.setdefault(key, i)
            if other != i:
                join.union(other, i)
            if len(lookup) > 1500000:
                raise MemoryError('Fixed 1.5 million exact-core limit; no eligibility result emitted.')
    return lookup, {key[:12] for key in lookup}


def matches(sequence, lookup, prefixes):
    reverse = rc(sequence)
    n = len(sequence)
    for j in range(n - 144):
        k = n - j - 145
        if sequence[j:j + 12] not in prefixes and reverse[k:k + 12] not in prefixes:
            continue
        hit = lookup.get(min(sequence[j:j + 145], reverse[k:k + 145]))
        if hit is not None:
            yield hit


def write_table(path, records):
    with gzip.open(path, 'wt') as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main():
    start = time.monotonic()
    prior = json.loads((OUT / 'additional_activity_spatial_and_reference_receipt.json').read_text())
    table = OUT / 'additional_activity_spatial_and_reference_screen.csv.gz'
    assert sha(table) == prior['screen_table_sha256']
    screened = list(rows(table))
    original = json.loads((OUT / 'additional_neuronal_activity_inventory_receipt.json').read_text())
    metadata = OUT / Path(original['sources'][0]['path']).name
    assert sha(metadata) == prior['metadata_sha256']
    # Keep only identity/location columns; do not load activity or P-value columns.
    fields = ['name', 'sequence', 'ref', 'chr', 'start', 'end', 'category']
    source = [{k: r[k] for k in fields} for r in rows(metadata, '\t')]
    by_name = {r['name']: i for i, r in enumerate(source)}
    assert len(by_name) == len(source) == prior['source_rows_before_measurement_filter']
    forward = LiftOver(str(ROOT / 'work/hg19ToHg38.over.chain.gz'))
    backward = LiftOver(str(ROOT / 'work/hg38ToHg19.over.chain.gz'))
    initial, mapped = independent_initial_components(source, forward, backward)
    assert len(set(initial)) == prior['all_source_spatial_and_full_DNA_components']
    for row in screened:
        i = by_name[row['name']]
        assert initial[i] == int(row['source_component'])
        coordinate = ':'.join(map(str, mapped[i])) if i in mapped else ''
        assert coordinate == row['mapped_GRCh38_interval']
    candidate_ids = {by_name[r['name']] for r in screened if r['passes_spatial_and_reference_screen'] == 'True'}
    unique_source = {}
    oriented_counts = collections.Counter()
    for i, row in enumerate(source):
        s = row['sequence'].upper()
        oriented_counts[s] += 1
        if not set(s) - set('ACGT'):
            unique_source.setdefault(min(s, rc(s)), i)
    join = Join(initial)
    iterations = []
    for iteration in range(1, 21):
        before = join.merges
        lookup, prefixes = build_lookup(unique_source, join, candidate_ids)
        nkeys = len(lookup)
        for sequence, i in unique_source.items():
            for other in matches(sequence, lookup, prefixes):
                join.union(i, other)
        iterations.append({'iteration': iteration, 'exact_canonical_145bp_cores': nkeys,
                           'additional_component_joins': join.merges - before})
        del lookup, prefixes
        print('EXACT CORE SOURCE CLOSURE', json.dumps(iterations[-1]), flush=True)
        if join.merges == before:
            break
    else:
        raise RuntimeError('No fixed point within 20 source-closure passes; no eligibility result emitted.')
    # Independently repeat spatial queries using a range-count formulation.
    history_points = collections.defaultdict(set)
    for entry in prior['historical_sources']:
        assert sha(entry['path']) == entry['sha256']
        for r in rows(entry['path']):
            history_points[r['chromosome']].add(int(r['position_hg38']) - 1)
    reserved_path = CHANNEL / 'GSE87711_design.txt.gz'
    assert sha(reserved_path) == prior['reserved_original_design_sha256']
    reserved_points = collections.defaultdict(set)
    reserved_sequences = set()
    with gzip.open(reserved_path, 'rt') as f:
        variants = set()
        for identifier, sequence in csv.reader(f, delimiter='\t'):
            variants.add(identifier.split('_')[0])
            reserved_sequences.add(min(sequence, rc(sequence)))
    for variant in variants:
        c, p, ref, alt = variant.split(':')
        for chrom, q, strand, score in forward.convert_coordinate('chr' + c.removeprefix('chr'), int(p) - 1):
            reserved_points[chrom].add(int(q))
    history_points = {c: sorted(v) for c, v in history_points.items()}
    reserved_points = {c: sorted(v) for c, v in reserved_points.items()}
    initial_history, initial_reserved = set(), set()
    for i, (c, a, b) in mapped.items():
        for points, target, margin in [(history_points, initial_history, 2001), (reserved_points, initial_reserved, 2146)]:
            values = points.get(c, [])
            if bisect.bisect_right(values, b + margin) - bisect.bisect_left(values, a - margin) > 0:
                target.add(initial[i])
    for r in screened:
        assert ('whole_source_component_near_any_historical_region' in r['screen_exclusions']) == (int(r['source_component']) in initial_history)
        assert ('whole_source_component_near_entire_reserved_study' in r['screen_exclusions']) == (int(r['source_component']) in initial_reserved)
    bad_history = {join.find(i) for i in initial_history}
    bad_reserved = {join.find(i) for i in initial_reserved}
    lookup, prefixes = build_lookup(unique_source, join, candidate_ids, bad_history | bad_reserved)
    literal_history, literal_reserved = set(), set()
    seen_history = set()
    history_profiles = []
    for entry in prior['historical_sources']:
        newly_scanned, hits = 0, set()
        for row in rows(entry['path']):
            for name in ['reference', 'alternate']:
                sequence = row[name].upper()
                if set(sequence) - set('ACGT'):
                    continue
                sequence = min(sequence, rc(sequence))
                digest = hashlib.sha256(sequence.encode()).digest()
                if digest in seen_history:
                    continue
                seen_history.add(digest)
                newly_scanned += 1
                hits.update(join.find(i) for i in matches(sequence, lookup, prefixes))
        literal_history.update(hits)
        history_profiles.append({'path': entry['path'], 'new_unique_full_sequences_scanned': newly_scanned,
                                 'candidate_components_with_exact_145bp_history_matches': len(hits)})
        print('EXACT CORE HISTORY SCREEN', json.dumps(history_profiles[-1]), flush=True)
    for sequence in reserved_sequences:
        literal_reserved.update(join.find(i) for i in matches(sequence, lookup, prefixes))
    del lookup, prefixes, seen_history
    records = []
    exclusions = collections.Counter()
    for row in screened:
        if row['passes_spatial_and_reference_screen'] != 'True':
            continue
        i = by_name[row['name']]
        root = join.find(i)
        reasons = []
        for bad, reason in [(bad_history, 'shared_core_linked_component_near_history'),
                            (bad_reserved, 'shared_core_linked_component_near_reserved_study'),
                            (literal_history, 'linked_component_shares_exact_145bp_with_history'),
                            (literal_reserved, 'linked_component_shares_exact_145bp_with_reserved_study')]:
            if root in bad:
                reasons.append(reason)
        if oriented_counts[source[i]['sequence'].upper()] != 1:
            reasons.append('same_oriented_DNA_has_multiple_original_construct_identifiers')
        exclusions.update(reasons)
        records.append({**row, 'source_component_after_shared_145bp_closure': root,
                        'passes_shared_core_and_duplicate_screen': not reasons,
                        'additional_exclusions': '|'.join(reasons),
                        'remaining': 'Independent no-prefix exact-core replay and a prespecified future activity experiment; current fit unchanged.'})
    target = OUT / 'additional_activity_closed_core_candidate_screen.csv.gz'
    write_table(target, records)
    full = OUT / 'all_neuronal_source_identity_components.csv.gz'
    write_table(full, [{'name': r['name'], 'source_index': i, 'initial_component': initial[i],
                       'component_after_shared_145bp_closure': join.find(i),
                       'mapped_GRCh38_interval': ':'.join(map(str, mapped[i])) if i in mapped else ''}
                      for i, r in enumerate(source)])
    kept = [r for r in records if r['passes_shared_core_and_duplicate_screen']]
    result = {'passed_exact_core_source_and_history_screen': True,
        'independent_initial_component_graph_and_spatial_queries_reproduced': True,
        'independent_replay_new_literal_rows': len(screened), 'all_source_rows': len(source),
        'exact_core_policy': 'Full 145bp string equality in either orientation. Iterate component expansion to a fixed point across the entire source, then protect all historical DNA and all original reserved-study constructs. Prefixes only accelerate exact string lookup.',
        'source_closure_iterations': iterations, 'source_closure_additional_joins': join.merges,
        'history_literal_profiles': history_profiles, 'reserved_full_construct_sequences_scanned': len(reserved_sequences),
        'source_original_build_and_reference_screen_sha256': sha(table),
        'metadata_sha256': sha(metadata), 'reserved_original_design_sha256': sha(reserved_path),
        'exclusion_counts_overlap': dict(exclusions), 'potential_rows_before': len(records),
        'potential_rows_after': len(kept),
        'potential_distinct_canonical_full_sequences_after': len({r['physical_sequence_sha256_both_orientations'] for r in kept}),
        'potential_components_after': len({r['source_component_after_shared_145bp_closure'] for r in kept}),
        'duplicate_policy': 'Exclude all original identifiers with identical oriented full DNA, regardless of coefficient values or availability. Opposite assay orientations remain distinct measurements; they share an exclusion component and are never averaged here.',
        'source_coefficients_or_P_used_for_exclusions': False, 'reserved_measurement_files_read': False,
        'models_fitted': False, 'running_comparison_changed': False, 'training_authorized': False,
        'candidate_screen_sha256': sha(target), 'all_source_components_sha256': sha(full),
        'remaining': ['Independent no-prefix shared-core replay on retained components against the entire source, historical DNA and entire reserved study',
                      'Separately specified activity experiment if current comparison warrants it; no new independent-study claim from this existing neuron study'],
        'code_sha256': sha(Path(__file__)), 'elapsed_seconds': time.monotonic() - start}
    (OUT / 'additional_activity_closed_core_receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()

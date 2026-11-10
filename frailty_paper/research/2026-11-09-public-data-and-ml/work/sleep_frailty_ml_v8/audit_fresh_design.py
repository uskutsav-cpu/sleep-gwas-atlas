"""Reserve source-measured physical DNA; never calculate model performance."""
import collections, re, bisect, importlib.util
from v8_common import *
sys.path.append(str(ROOT/'work/stage5_env/lib/python3.13/site-packages'))
from pyliftover import LiftOver
# Read-only worksheet helper has no dependency on the old model runtime.
spec = importlib.util.spec_from_file_location('readonly_source_xlsx', ROOT/'work/sleep_frailty_ml_v7/extract_neuroimmune_sources.py')
reader = importlib.util.module_from_spec(spec); spec.loader.exec_module(reader)
spec = importlib.util.spec_from_file_location('official_twobit_reader', ROOT/'work/sleep_frailty_ml_v7/reference.py')
# reference.py imports "common"; do not call its receipt-writing get_reference.
sys.path.append(str(ROOT/'work/sleep_frailty_ml_v7'))
refmodule = importlib.util.module_from_spec(spec); spec.loader.exec_module(refmodule)


def main():
    assert not (OUT/'full_model_search_started.json').exists()
    assert_previous_intact()
    source = V7/'sources/S12C1A15T_PMC11839127_media-3_ENIGMA_MPRA_results.xlsx'
    tables = reader.read_tables(source, ['Legends', 'a. Regulatory activity (Veh)', 'b. Allelic effect (Veh)'])
    legends = tables['Legends']['rows']
    assert any('hg19' in str(row) for row in legends)
    raw = tables['a. Regulatory activity (Veh)']['rows']
    header = raw[0]; metadata = collections.defaultdict(set)
    for row in raw[1:]:
        d = dict(zip(header, row))
        if d['SNP'].startswith('rs'):
            metadata[d['SNP']].add(('chr'+d['chr'].removeprefix('chr'), int(float(d['pos']))))
    b = tables['b. Allelic effect (Veh)']['rows']
    # Store only identities and availability, not coefficient, p or q values.
    effect_ids = [row[0] for row in b[1:]]
    assert len(effect_ids) == len(set(effect_ids))
    effect_metadata = {row[0]: {'source_ref_label': row[1], 'source_alt_label': row[2]} for row in b[1:]}
    effect_status = dict(collections.Counter(row[-1] for row in b[1:]))
    del tables, raw, b
    physical = {}; repeats = collections.Counter()
    with gzip.open(SOURCES/'GSE289858_variant_barcode_statistics.txt.gz', 'rt') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            name, seq = row['name'], row['variant'].upper()
            if name in physical:
                assert physical[name] == seq
                repeats[name] += 1
            else:
                physical[name] = seq
    fasta_url = 'https://api.bitbucket.org/2.0/repositories/steinlabunc/enigma-mpra/src/main/data/release/ENIGMA_LD_noprimers.fasta'
    fasta_receipt = fetch(fasta_url, SOURCES/'author_ENIGMA_LD_noprimers.fasta', 10_000_000)
    author = {}; name = None
    for line in (SOURCES/'author_ENIGMA_LD_noprimers.fasta').read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].strip(); assert name not in author; author[name] = ''
        else:
            author[name] += line.strip().upper()
    shared = set(author)&set(physical)
    # GEO's barcode mapping is stored in the reverse-complement sequencing
    # orientation. All shared 150-bp DNA strings match the author design after
    # this explicit transform. Preserve the author's construct orientation;
    # do not silently assume reporter effects are orientation-invariant.
    assert len(shared) > 15000 and all(author[k] == rc(physical[k]) for k in shared)
    original_GEO_count = len(physical)
    physical = {k: author[k] for k in shared}
    genome_path = V6/'sources/official_hg38.2bit'
    ref_audit = json.loads((V6/'outputs/additional_reference_audit.json').read_text())
    assert sha(genome_path) == ref_audit['sha256']
    genome = refmodule.TwoBit(genome_path)
    lift = LiftOver(str(ROOT/'work/hg19ToHg38.over.chain.gz'))
    back = LiftOver(str(ROOT/'work/hg38ToHg19.over.chain.gz'))
    groups = collections.defaultdict(dict)
    excluded = collections.Counter()
    for name, seq in physical.items():
        m = re.search(r'_SNP_(rs\d+)_.*_A([12])(?:_|$)', name)
        if not m:
            excluded['control_or_name_without_unambiguous_rsid_and_A1_A2'] += 1; continue
        normalized = re.sub(r'_A[12](?=_|$)', '_ALLELE', name)
        assert m[2] not in groups[normalized]
        groups[normalized][m[2]] = (name, seq, m[1])
    verified = []; by_rsid = collections.defaultdict(list)
    for identity, group in groups.items():
        if set(group) != {'1', '2'}:
            excluded['unpaired_design'] += 1; continue
        an, a, rsid = group['1']; bn, bb, rsid2 = group['2']
        assert rsid == rsid2
        if len(a) != len(bb) or len(a) != 150 or set(a+bb)-set('ACGT'):
            excluded['not_unambiguous150bp_pair'] += 1; continue
        diff = [i for i, (u, v) in enumerate(zip(a, bb)) if u != v]
        if len(diff) != 1:
            excluded['not_actual_single_SNV'] += 1; continue
        coords = metadata.get(rsid, set())
        if len(coords) != 1:
            excluded['missing_or_ambiguous_original_hg19_coordinate'] += 1; continue
        chrom, pos = next(iter(coords)); p = pos-1
        maps = lift.convert_coordinate(chrom, p)
        if len(maps) != 1:
            excluded['nonunique_liftover'] += 1; continue
        c, q, strand, _ = maps[0]
        reverse = back.convert_coordinate(c, q)
        if len(reverse) != 1 or reverse[0][:3] != (chrom, p, strand):
            excluded['nonreciprocal_liftover'] += 1; continue
        j = diff[0]; matches = []
        for orientation in ['+', '-']:
            at = j if orientation == '+' else len(a)-1-j
            start = q-at
            if start < 0:
                continue
            truth = genome.sequence(c, start, start+len(a))
            for reference_is_A1, candidate in [(True, a), (False, bb)]:
                if truth == (candidate if orientation == '+' else rc(candidate)):
                    matches.append((orientation, reference_is_A1, start))
        if len(matches) != 1:
            excluded['full_physical_reference_genome_not_unique_match'] += 1; continue
        orientation, a1ref, start = matches[0]
        ref, alt = (a, bb) if a1ref else (bb, a)
        ra, aa = (ref[j], alt[j]) if orientation == '+' else (rc(ref[j]), rc(alt[j]))
        row = {'rsid': rsid, 'variant_GRCh38': f'{c}:{q+1}:{ra}:{aa}', 'chromosome': c, 'position_hg38': q+1, 'original_hg19_chromosome': chrom, 'original_hg19_position1': pos, 'liftover_strand': strand, 'construct_strand': orientation, 'reference_name': an if a1ref else bn, 'alternate_name': bn if a1ref else an, 'reference': ref, 'alternate': alt, 'changed_position': j, 'full_reference_matches_hg38': True, 'source_baseline_allelic_coefficient_available': rsid in effect_metadata}
        by_rsid[rsid].append(row)
    for rsid, rows in by_rsid.items():
        identities = {(r['reference'], r['alternate'], r['variant_GRCh38']) for r in rows}
        if len(identities) != 1:
            excluded['multiple_distinct_physical_contexts_per_rsid'] += len(rows); continue
        verified.append(rows[0])
        excluded['identical_physical_pairs_with_different_locus_names_collapsed'] += len(rows)-1
    assert len(verified) > 5000, (len(verified), excluded)
    verified.sort(key=lambda r: (r['chromosome'], r['position_hg38'], r['variant_GRCh38']))
    write_csv(DATA/'matoba_full_verified_physical_design.csv', verified)
    # Exclude the whole fresh 2-kb component if it touches either revealed V7
    # outside experiment. Source rows are NOT filtered by allelic significance.
    known_positions = collections.defaultdict(list); known_cores = set()
    for path in [V7/'reserved_physical_design.csv', V7/'neuroimmune_physical_design_verified.csv']:
        for r in csv_rows(path):
            known_positions[r['chromosome']].append(int(r['position_hg38']))
            a, bb = r['reference'], r['alternate']; j = int(r['changed_position'])
            lo = max(0, j-149); hi = min(j, len(a)-150)
            for s in [a, bb]:
                for k in range(lo, hi+1):
                    core = s[k:k+150]; known_cores.add(min(core, rc(core)))
    for c in known_positions:
        known_positions[c].sort()
    groups = []; current = []
    for row in verified:
        if current and (row['chromosome'] != current[-1]['chromosome'] or row['position_hg38']-current[-1]['position_hg38'] > 2000):
            groups.append(current); current = []
        current.append(row)
    if current:
        groups.append(current)
    allowed = []
    for rows in groups:
        bad = False
        for r in rows:
            ps = known_positions[r['chromosome']]; j = bisect.bisect_left(ps, r['position_hg38'])
            if any(0 <= k < len(ps) and abs(ps[k]-r['position_hg38']) <= 2000 for k in [j-1, j]):
                bad = True
            if any(min(s, rc(s)) in known_cores for s in [r['reference'], r['alternate']]):
                bad = True
        if bad:
            excluded['new_component_overlaps_revealed_prior_outside_neighborhood_or_core'] += len(rows)
        else:
            allowed.extend(rows)
    # Audit exact native 150-bp contexts against ALL public pretraining rows
    # because the example checkpoint's original fold is not specified.
    candidate = collections.defaultdict(set)
    for r in allowed:
        for s in [r['reference'], r['alternate']]:
            candidate[min(s, rc(s))].add(r['rsid'])
    prefixes = {s[:16] for s in candidate}; pretraining_seen = set(); npre = 0
    with (SOURCES/'legnet_datasets_original_WTC11.tsv').open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            s = r['seq']; npre += 1
            for u in [s, rc(s)]:
                for k in range(len(u)-149):
                    if u[k:k+16] in prefixes:
                        core = u[k:k+150]
                        if core in candidate:
                            pretraining_seen.update(candidate[core])
    allowed = [r for r in allowed if r['rsid'] not in pretraining_seen]
    excluded['exact150bp_pretraining_overlap'] = len(pretraining_seen)
    write_csv(DATA/'reserved_matoba_physical_design.csv', allowed)
    # Future V8 training removes ENTIRE V7 groups touching these designs.
    points = collections.defaultdict(list); test_cores = set()
    for r in allowed:
        points[r['chromosome']].append(r['position_hg38'])
        for s in [r['reference'], r['alternate']]:
            test_cores.add(min(s, rc(s)))
    for c in points:
        points[c].sort()
    excluded_training = set()
    for r in csv_rows(V7/'sequence_bank.csv.gz'):
        ps = points[r['chromosome']]; p = int(r['position_hg38']); j = bisect.bisect_left(ps, p)
        bad = any(0 <= k < len(ps) and abs(ps[k]-p) <= 2000 for k in [j-1, j])
        if not bad:
            a, bb = r['reference'], r['alternate']; at = int(r['changed_position'])
            for s in [a, bb]:
                for k in range(max(0, at-149), min(at, len(a)-150)+1):
                    core = s[k:k+150]
                    if min(core, rc(core)) in test_cores:
                        bad = True; break
                if bad:
                    break
        if bad:
            excluded_training.add(r['locus_group'])
    write_json(OUT/'V7_development_groups_excluded_for_fresh_test.json', sorted(excluded_training))
    # Keep source identities without effect values for the locked evaluator.
    write_json(OUT/'reserved_matoba_effect_identity_metadata.json', {r['rsid']: effect_metadata[r['rsid']] for r in allowed if r['rsid'] in effect_metadata})
    receipt = {
        'original_source_workbook_sha256': sha(source), 'source_workbook_version': 'Archived original preprint supplement; latest journal-supplement equivalence not yet established',
        'source_coordinate_build_verified_from_legend': 'hg19/GRCh37, 1-based',
        'published_baseline_original_coefficients': len(effect_ids), 'published_detection_flags_unfiltered': effect_status,
        'baseline_source_scope': 'Original MPRALM coefficients only for author-defined active elements, including both detected and undetected allelic effects. Not a representative all-variant or patient sample.',
        'author_fasta_receipt': fasta_receipt, 'GEO_and_author_fasta_exact_reverse_complement_shared_constructs': len(shared),
        'orientation_policy': 'Preserve author design orientation after exact reverse-complement reconciliation with the barcode table. This is source reconciliation, not a claim of assay orientation invariance.',
        'physical_name_duplicates': dict(repeats), 'unique_GEO_constructs': original_GEO_count,
        'GEO_constructs_without_matching_author_design_excluded': original_GEO_count-len(shared),
        'fully_genome_verified_unique_physical_pairs': len(verified), 'reserved_eligible_physical_pairs': len(allowed),
        'reserved_pairs_with_source_original_baseline_coefficient': sum(r['rsid'] in effect_metadata for r in allowed),
        'all_pretraining_source_rows_scanned': npre, 'pretraining_exact150bp_overlap_removed': len(pretraining_seen),
        'pretraining_full_genomic_neighborhood_exclusion_established': False,
        'V7_development_whole_groups_to_exclude_before_V8_fit': len(excluded_training),
        'exclusions': dict(excluded), 'original_fitted_coefficients_or_pvalues_exported_here': False,
        'fresh_model_performance_calculated': False, 'single_donor_scope': 'Primary human neural progenitors, male donor 88, repeated cultures; donor independence from prior work not established.',
        'unresolved': ['Journal vs preprint supplement equivalence', 'Original reported signed allele convention, beyond magnitude and genome-pair verification', 'Source raw barcode/model reconstruction not independently performed', 'One donor and active-element selection limit generalizability'],
        'outputs_sha256': {n: sha(DATA/n) for n in ['matoba_full_verified_physical_design.csv', 'reserved_matoba_physical_design.csv']},
        'reader_source_sha256': sha(ROOT/'work/sleep_frailty_ml_v7/reference.py'),
    }
    write_json(OUT/'fresh_matoba_physical_design_audit.json', receipt)
    assert_previous_intact()
    print('FRESH STUDY PHYSICAL-DNA AUDIT COMPLETE', json.dumps(receipt), flush=True)


if __name__ == '__main__':
    main()

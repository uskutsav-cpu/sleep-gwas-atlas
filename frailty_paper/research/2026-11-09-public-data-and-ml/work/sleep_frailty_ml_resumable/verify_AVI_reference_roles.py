"""Resolve source allele roles using complete actual DNA and human references."""
import ast
import collections
import concurrent.futures
import csv
import gzip
import hashlib
import json
from pathlib import Path
import ssl
import struct
import sys
import time
import urllib.parse
import urllib.request
import certifi
from audit_permitted_AVI_alignment import ROOT, BASE, DATA, OUT, sha

sys.path.append(str(ROOT / 'work/stage5_env/lib/python3.13/site-packages'))
from pyliftover import LiftOver


def rc(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def main():
    start = time.monotonic()
    audit_path = OUT / 'permitted_AVI_alignment_receipt.json'
    audit = json.loads(audit_path.read_text())
    cases = {r['bank_id']: r for r in audit['failures'] if r['reason'] == 'reference_allele_disagrees_with_official_GRCh38'}
    assert len(cases) == 601 and all(r['genome_base'] == r['variant'].split(':')[3] for r in cases.values())
    reader_path = ROOT / 'work/sleep_frailty_ml_v7/reference.py'
    assert sha(reader_path) == audit['reference_reader_sha256']
    tree = ast.parse(reader_path.read_text())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'TwoBit')
    ns = {'struct': struct}
    exec(compile(ast.Module(body=[cls], type_ignores=[]), str(reader_path), 'exec'), ns)
    genome_path = BASE / 'sleep_frailty_ml_v6/sources/official_hg38.2bit'
    assert sha(genome_path) == audit['genome_sha256']
    genome = ns['TwoBit'](genome_path)
    metadata = json.loads((ROOT / 'outputs/sleep_frailty_ml_resumable/development/frozen_reader_none_42_started.json').read_text())
    raw = {}
    for entry in metadata['physical_source_receipt']['sources']:
        assert sha(Path(entry['path'])) == entry['sha256']
        with gzip.open(entry['path'], 'rt') as f:
            for row in csv.DictReader(f):
                i = int(row['bank_id'])
                if i in cases:
                    assert i not in raw
                    raw[i] = row
    assert raw.keys() == cases.keys()
    chain_receipt = json.loads((ROOT / 'outputs/sleep_frailty_ml_neuronal_activity/additional_activity_spatial_and_reference_receipt.json').read_text())
    fpath, bpath = ROOT / 'work/hg19ToHg38.over.chain.gz', ROOT / 'work/hg38ToHg19.over.chain.gz'
    for p in [fpath, bpath]:
        assert sha(p) == chain_receipt['chain_sha256'][p.name]
    forward, backward = LiftOver(str(fpath)), LiftOver(str(bpath))
    results, groups = {}, collections.defaultdict(list)
    for i, row in raw.items():
        chrom, pos, a, b = row['variant_GRCh38'].split(':')
        p = int(pos) - 1
        actual_ref, actual_alt, j = row['reference'], row['alternate'], int(row['changed_position'])
        plus = actual_ref[j] + actual_alt[j] == a + b
        assert plus or rc(actual_ref[j]) + rc(actual_alt[j]) == a + b
        length = len(actual_ref)
        begin = p - (j if plus else length - 1 - j)
        truth = genome.sequence(chrom, begin, begin + length)
        rr, aa = (actual_ref, actual_alt) if plus else (rc(actual_ref), rc(actual_alt))
        assert truth[p - begin] == b
        result = {'bank_id': i, 'original_assay_variant_tag': row['variant_GRCh38'],
            'canonical_genomic_variant': f'{chrom}:{pos}:{b}:{a}',
            'GRCh38_reference_is_source_assay_allele': 'alternate',
            'full_construct_reference_verified': False, 'full_reference_match_assembly': '',
            'matched_assay_allele': '', 'source_assay_orientation_to_GRCh38': '+' if plus else '-',
            'GRCh38_span_start_0_based': begin, 'source_construct_length': length,
            'GRCh38_background_mismatches_for_assay_alternate': sum(x != y for x, y in zip(aa, truth)),
            'GRCh37_chromosome': '', 'GRCh37_variant_position_0_based': '',
            'GRCh37_span_start_0_based': '', 'GRCh37_source_assay_orientation': '',
            'reciprocal_unique_variant_build_mapping': False, 'reference_payload_sha256': '',
            'unresolved_reason': '', 'source_DNA_or_signed_targets_changed': False}
        if aa == truth:
            result.update(full_construct_reference_verified=True, full_reference_match_assembly='GRCh38', matched_assay_allele='alternate')
        else:
            maps = backward.convert_coordinate(chrom, p)
            if len(maps) != 1:
                result['unresolved_reason'] = 'nonunique_or_missing_GRCh38_to_GRCh37_mapping'
            else:
                c19, p19, strand, score = maps[0]
                back = forward.convert_coordinate(c19, p19)
                if len(back) != 1 or back[0][:3] != (chrom, p, strand):
                    result['unresolved_reason'] = 'nonreciprocal_variant_build_mapping'
                else:
                    old_plus = plus == (strand == '+')
                    begin19 = int(p19) - (j if old_plus else length - 1 - j)
                    result.update(GRCh37_chromosome=c19, GRCh37_variant_position_0_based=int(p19),
                        GRCh37_span_start_0_based=begin19, GRCh37_source_assay_orientation='+' if old_plus else '-',
                        reciprocal_unique_variant_build_mapping=True)
                    if begin19 < 0:
                        result['unresolved_reason'] = 'negative_GRCh37_construct_start'
                    else:
                        groups[(c19, begin19 // 16384)].append(i)
        results[i] = result
    cache = OUT / 'GRCh37_reference_slices'
    cache.mkdir(exist_ok=True)
    context = ssl.create_default_context(cafile=certifi.where())
    def query(item):
        (chrom, window), ids = item
        lo = min(results[i]['GRCh37_span_start_0_based'] for i in ids)
        hi = max(results[i]['GRCh37_span_start_0_based'] + results[i]['source_construct_length'] for i in ids)
        args = urllib.parse.urlencode({'genome': 'hg19', 'chrom': chrom, 'start': lo, 'end': hi})
        url = 'https://api.genome.ucsc.edu/getData/sequence?' + args
        path = cache / (hashlib.sha256(url.encode()).hexdigest()[:24] + '.json')
        if path.exists():
            content = path.read_bytes()
        else:
            for attempt in range(3):
                try:
                    with urllib.request.urlopen(url, context=context, timeout=20) as response:
                        content = response.read(100001)
                    assert len(content) <= 100000
                    data = json.loads(content)
                    assert data['genome'] == 'hg19' and data['chrom'] == chrom and data['start'] == lo and data['end'] == hi
                    assert len(data['dna']) == hi - lo
                    tmp = path.with_suffix('.pending')
                    tmp.write_bytes(content)
                    tmp.replace(path)
                    break
                except Exception:
                    if attempt == 2:
                        raise
                    time.sleep(1 + attempt)
        data = json.loads(content)
        assert data['genome'] == 'hg19' and data['chrom'] == chrom and data['start'] == lo and data['end'] == hi
        dna = data['dna'].upper()
        assert len(dna) == hi - lo and not set(dna) - set('ACGTN')
        return ids, lo, dna, {'url': url, 'path': str(path.relative_to(ROOT)), 'sha256': sha(path)}
    failures, payloads, completed = [], [], 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        jobs = {pool.submit(query, item): item[0] for item in groups.items()}
        for job in concurrent.futures.as_completed(jobs):
            completed += 1
            try:
                ids, lo, dna, payload = job.result()
                payloads.append(payload)
                for i in ids:
                    result = results[i]
                    off = result['GRCh37_span_start_0_based'] - lo
                    truth = dna[off:off + result['source_construct_length']]
                    r = raw[i]
                    plus = result['GRCh37_source_assay_orientation'] == '+'
                    a, b = (r['reference'], r['alternate']) if plus else (rc(r['reference']), rc(r['alternate']))
                    matched = [label for label, s in [('reference', a), ('alternate', b)] if s == truth]
                    result['reference_payload_sha256'] = payload['sha256']
                    if len(matched) == 1:
                        result.update(full_construct_reference_verified=True, full_reference_match_assembly='GRCh37', matched_assay_allele=matched[0])
                    else:
                        result['unresolved_reason'] = 'neither_complete_assay_allele_matches_GRCh37_or_GRCh38_reference'
            except Exception as error:
                failures.append({'group': jobs[job], 'error': repr(error)})
            if completed % 10 == 0 or completed == len(groups):
                print('PUBLIC GRCh37 REFERENCE WINDOWS', completed, '/', len(groups), 'failed', len(failures), flush=True)
    target = OUT / 'canonical_reference_role_resolution.csv'
    with target.open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(next(iter(results.values()))))
        writer.writeheader()
        writer.writerows(results[i] for i in sorted(results))
    counts = collections.Counter(r['full_reference_match_assembly'] or 'unresolved' for r in results.values())
    receipt = {'reference_resolution_complete_without_network_errors': not failures,
        'source_reversed_reference_constructs': len(results), 'full_construct_match_counts': dict(counts),
        'unresolved_reasons': dict(collections.Counter(r['unresolved_reason'] for r in results.values() if not r['full_construct_reference_verified'])),
        'network_failures': failures, 'source_alignment_audit_sha256': sha(audit_path),
        'genome_sha256': audit['genome_sha256'], 'public_GRCh37_reference_payloads': payloads,
        'chain_sha256': {p.name: sha(p) for p in [fpath, bpath]},
        'policy': 'Use a canonical GRCh38 reference-to-alternate score key only after a whole actual construct matches GRCh38 or matches GRCh37 at a unique reciprocal mapped variant position. Keep unresolved backgrounds unavailable; never edit source DNA or invert experimental signed targets.',
        'original_genetic_discovery_effect_allele_conventions_established': False,
        'original_source_DNA_targets_or_models_changed': False, 'reserved_outcomes_read': False,
        'resolution_table_sha256': sha(target), 'code_sha256': sha(Path(__file__)),
        'elapsed_seconds': time.monotonic() - start}
    (OUT / 'canonical_reference_role_resolution_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ['source_reversed_reference_constructs', 'full_construct_match_counts', 'unresolved_reasons', 'network_failures', 'elapsed_seconds']}), flush=True)
    assert not failures, failures


if __name__ == '__main__':
    main()

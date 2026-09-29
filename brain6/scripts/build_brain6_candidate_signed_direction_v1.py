#!/usr/bin/env python3
"""Look up the frozen candidate lead signs in receipt-bound five-pair PLACO output."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CAND = ROOT / 'brain6/results/loci/placo_five_track_candidate_loci_v1/candidate_loci.tsv'
MASTER = ROOT / 'brain6/results/track_b_admission_v1/placo_master_five_track_v1.tsv'
OUT = ROOT / 'brain6/results/brain6_candidate_signed_direction_v1'


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def write(path: Path, rows: list[dict]) -> None:
    with path.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, delimiter='\t', fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    if OUT.exists():
        raise FileExistsError('Signed-direction output exists; refuse overwrite')
    candidates = read(CAND)
    provenance = json.loads((CAND.parent / 'provenance.json').read_text())
    if len(candidates) != 25 or sha(CAND) != provenance['candidate_loci_sha256']:
        raise ValueError('Frozen 25-candidate file changed')
    master = read(MASTER)
    by_pair = {r['pair_id']: r for r in master}
    if len(by_pair) != 5:
        raise ValueError('Expected five source pair outputs')
    targets = defaultdict(set)
    for candidate in candidates:
        for lead in candidate['lead_variants'].split(';'):
            targets[candidate['pair_id']].add(lead)
    found = {}
    source_sha = {}
    for pair, leads in targets.items():
        meta = by_pair[pair]
        if not meta['output_path'].startswith('brain6/results/placo/'):
            raise ValueError('Unexpected PLACO output path')
        path = ROOT / meta['output_path']
        digest = sha(path)
        if digest != meta['output_sha256'] or digest != provenance['pair_output_sha256'][pair]:
            raise ValueError(f'Frozen PLACO pair output changed: {pair}')
        source_sha[str(path.relative_to(ROOT))] = digest
        with gzip.open(path, 'rt', newline='') as stream:
            reader = csv.DictReader(stream, delimiter='\t')
            for row in reader:
                key = (pair, row['SNP'])
                if row['SNP'] in leads:
                    if key in found:
                        raise ValueError(f'Duplicate source lead row: {key}')
                    found[key] = row
    detailed = []
    summary = []
    for candidate in candidates:
        pair = candidate['pair_id']
        relations = []
        products = []
        for lead in candidate['lead_variants'].split(';'):
            source = found.get((pair, lead))
            if source is None or source['status'] != 'TESTED':
                raise ValueError(f'Missing tested frozen PLACO lead: {pair}/{lead}')
            if source['CHR'] != candidate['CHR'] or not int(candidate['START']) <= int(source['BP']) <= int(candidate['STOP']):
                raise ValueError(f'Lead coordinate mismatch: {pair}/{lead}')
            z1, z2 = float(source['Z1']), float(source['Z2'])
            relation = 'CONCORDANT' if z1 * z2 > 0 else 'DISCORDANT' if z1 * z2 < 0 else 'ZERO_SIGN'
            relations.append(relation)
            products.append(z1 * z2)
            detailed.append({'candidate_locus_id': candidate['locus_id'], 'pair_id': pair,
                             'lead_variant': lead, 'chromosome': source['CHR'], 'position_grch37': source['BP'],
                             'z_sleep': source['Z1'], 'z_disorder': source['Z2'],
                             'signed_z_product': z1 * z2, 'direction_relation': relation,
                             'placo_p': source['P_PLACO'], 'analysis_label': 'DESCRIPTIVE_SOURCE_Z_DIRECTION'})
        summary.append({'candidate_locus_id': candidate['locus_id'], 'pair_id': pair,
                        'lead_variants': candidate['lead_variants'], 'n_leads': len(relations),
                        'direction_relations': ';'.join(relations),
                        'all_leads_concordant': all(r == 'CONCORDANT' for r in relations),
                        'any_lead_discordant': any(r == 'DISCORDANT' for r in relations),
                        'signed_z_products': ';'.join(str(x) for x in products),
                        'interpretation': 'Source Z signs only; not local covariance or causal direction'})
    if len(summary) != 25 or len(detailed) != 27:
        raise ValueError('Lead coverage changed')
    OUT.mkdir(parents=True)
    detail_path = OUT / 'lead_directions.tsv'
    summary_path = OUT / 'candidate_directions_25.tsv'
    write(detail_path, detailed)
    write(summary_path, summary)
    receipt = {'schema_version': 1, 'candidate_count': len(summary), 'lead_count': len(detailed),
               'candidate_sha256': sha(CAND), 'master_sha256': sha(MASTER),
               'pair_source_sha256': source_sha,
               'lead_directions_sha256': sha(detail_path),
               'candidate_directions_sha256': sha(summary_path),
               'builder_sha256': sha(Path(__file__)),
               'interpretation': 'Descriptive signed Z relation from original five-pair PLACO outputs; no LAVA or source modification'}
    with (OUT / 'receipt.json').open('x') as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write('\n')
    print(json.dumps({'candidate_count': 25, 'lead_count': 27,
                      'concordant_candidates': sum(r['all_leads_concordant'] for r in summary)}))


if __name__ == '__main__':
    main()

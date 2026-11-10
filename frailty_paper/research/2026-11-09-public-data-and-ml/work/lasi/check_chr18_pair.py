#!/usr/bin/env python3
"""Check one prespecified candidate pair in the LASI-DAD public LD lookup."""
import csv, gzip, hashlib, json, os

source = 'work/lasi/chr18_ld.csv.gz'
expected_bytes = 442111058
if os.path.getsize(source) != expected_bytes:
    raise SystemExit('Incomplete download: refusing to interpret partial file')

target = {55206086, 55389957}  # provisional GRCh38 candidate coordinates
matches = []
digest = hashlib.sha256()
with open(source, 'rb') as fh:
    for block in iter(lambda: fh.read(1024 * 1024), b''):
        digest.update(block)
with gzip.open(source, 'rt', newline='') as fh:
    rows = csv.DictReader(fh)
    for row in rows:
        # Source positions can be written in scientific notation. The variant
        # identifiers retain the integer coordinates and allele contrast.
        pos1 = int(row['variant1'].split(':')[1])
        pos2 = int(row['variant2'].split(':')[1])
        if {pos1, pos2} == target:
            matches.append(row)

receipt = {
    'source_bytes': expected_bytes,
    'source_sha256': digest.hexdigest(),
    'source_url': 'https://github.com/zhengli09/LASI_DAD-LD-Reference-Panel',
    'candidates': ['rs1942262', 'rs11152363'],
    'provisional_grch38_positions': sorted(target),
    'matches': matches,
    'limitations': 'Pair lookup is restricted to r2 >= 0.2 and variants within 1 Mb; absence does not prove zero LD or presence of both variants. Verify alleles and genomic build. Haplotype correlation is not phenotype validation.'
}
with open('outputs/stage6_lasi_chr18_candidate_ld_receipt.json', 'w') as fh:
    json.dump(receipt, fh, indent=2)
print(json.dumps(receipt, indent=2))

#!/usr/bin/env python3
"""Count exact candidate alleles appearing in the filtered LASI-DAD LD pairs."""
import csv, gzip, hashlib, json, os

source = 'work/lasi/chr18_ld.csv.gz'
if os.path.getsize(source) != 442111058:
    raise SystemExit('Incomplete download')
targets = {'18:55206086': 'rs1942262', '18:55389957': 'rs11152363'}
counts = {v: {} for v in targets.values()}
with gzip.open(source, 'rt', newline='') as f:
    header = next(f)
    for line in f:
        if not any(k in line for k in targets):
            continue
        row = next(csv.reader([line]))
        for ix, afidx in ((5, 7), (6, 8)):
            variant = row[ix]
            key = ':'.join(variant.split(':')[:2])
            if key in targets:
                entries = counts[targets[key]]
                if variant not in entries:
                    entries[variant] = {'rows': 0, 'alternate_frequency': row[afidx]}
                entries[variant]['rows'] += 1
out = {'header': header.strip(), 'counts_by_candidate': counts,
       'limitations': 'This table contains only reported LD pairs r2 >= 0.2 within 1 Mb. Presence does not establish target cohort score eligibility; absence does not establish variant absence or r2=0. Coordinates provisionally GRCh38; inspect reference and source allele coding.'}
with open('outputs/stage6_lasi_chr18_marker_presence.json', 'w') as f:
    json.dump(out, f, indent=2)
print(json.dumps(out, indent=2))

"""Inventory direct project-candidate matches in the public PGS005221 file.

Matches are descriptive and do not test either score or infer biology.
"""
import csv, gzip, hashlib, json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
output, work = root / 'outputs', root / 'work'
sources = [('15', 'unvalidated_15_variant_score_prototype.csv', 'chr', 'position_GRCh37_inferred'),
           ('22', 'unvalidated_22_shared_direction_score.csv', 'chr_GRCh37_assumed', 'position_GRCh37_assumed')]
by_position, candidates = {}, {}
for label, name, chrom, pos in sources:
    with (output / name).open() as f:
        for row in csv.DictReader(f):
            key = (row[chrom], int(row[pos]))
            candidate = (label, row['variant'], row['counted_allele'], row['other_allele'])
            by_position.setdefault(key, []).append(candidate)
            candidates[(label, row['variant'])] = key
pgs_file = work / 'PGS005221.txt.gz'
matches, count = [], 0
with gzip.open(pgs_file, 'rt') as f:
    reader = csv.DictReader((line for line in f if not line.startswith('#')), delimiter='\t')
    for row in reader:
        count += 1
        key = (row['chr_name'], int(row['chr_position']))
        for label, rsid, counted, other in by_position.get(key, []):
            matches.append({'model': label, 'variant': rsid, 'chromosome': key[0], 'position37': key[1],
                            'pgs_rsid': row['rsID'], 'rsid_match': row['rsID'] == rsid,
                            'pgs_effect_allele': row['effect_allele'], 'pgs_other_allele': row['other_allele'],
                            'project_counted_allele': counted, 'project_other_allele': other,
                            'allele_pair_match': {counted, other} == {row['effect_allele'], row['other_allele']},
                            'pgs_effect_weight': row['effect_weight']})
if count != 943126:
    raise ValueError(f'Expected 943126 records, observed {count}')
fields = ['model', 'variant', 'chromosome', 'position37', 'pgs_rsid', 'rsid_match', 'pgs_effect_allele', 'pgs_other_allele', 'project_counted_allele', 'project_other_allele', 'allele_pair_match', 'pgs_effect_weight']
with (output / 'stage6_frailty_general_pgs_matches.csv').open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader(); writer.writerows(matches)
receipt = {'file': 'PGS005221.txt.gz', 'source': 'https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS005221/ScoringFiles/PGS005221.txt.gz',
           'sha256': hashlib.sha256(pgs_file.read_bytes()).hexdigest(), 'score_rows': count,
           'position_matches': {label: sum(row['model'] == label for row in matches) for label in ('15', '22')},
           'exact_position_and_pair_matches': {label: sum(row['model'] == label and row['allele_pair_match'] for row in matches) for label in ('15', '22')},
           'limitations': 'Position and unordered allele pairs are not sufficient to resolve multiallelic source contrasts; scores and outcome performance are not evaluated.'}
(output / 'stage6_frailty_general_pgs_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(receipt, indent=2))

"""Join two existing candidate screens by GRCh37 position and alleles.

The join is descriptive; the Lee estimates are an unvalidated reconstruction.
"""
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OTHER = HERE.parent / 'siraj_2026_mpra' / 'sknsh_results.tsv'
LEE = HERE / 'provisional_allelic_42.tsv'
OUT = HERE / 'cross_assay_candidate_comparison.tsv'

def read(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

lee, siraj = read(LEE), read(OTHER)
index = {}
for row in lee:
    chrompos, alleles = row['source_key'].split('_')
    key = chrompos + ':' + alleles
    if key in index:
        raise ValueError('duplicate Lee key: ' + key)
    index[key] = row

out = []
for s in siraj:
    l = index.get(s['variant37'])
    if not l:
        continue
    if l['rsid'] != s['rsid']:
        raise ValueError('rsID disagreement: ' + s['variant37'])
    out.append({
        'rsid': l['rsid'], 'variant37': s['variant37'], 'region': l['region'],
        'lee_listed_excluded': l['listed_excluded'],
        'lee_v124_provisional_logFC': l['v124_logFC'],
        'lee_v124_provisional_bh': l['v124_adj.P.Val'],
        'lee_v134_provisional_logFC': l['v134_logFC'],
        'lee_v134_provisional_bh': l['v134_adj.P.Val'],
        'siraj_library': s['library'], 'siraj_log2_skew': s['log2_skew'],
        'siraj_skew_se': s['skew_se'], 'siraj_source_fdr': s['source_fdr'],
        'siraj_A_activity_bonferroni_p': s['a_activity_bonferroni_p'],
        'siraj_B_activity_bonferroni_p': s['b_activity_bonferroni_p'],
    })
if len(out) != 4:
    raise ValueError('unexpected overlap; inspect before interpreting: ' + str(len(out)))
with OUT.open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(out[0]), delimiter='\t')
    writer.writeheader(); writer.writerows(out)
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
receipt = {'join': 'exact GRCh37 chr:position:REF:ALT and identical rsID',
           'rows': len(out), 'unique_variants': len({x['rsid'] for x in out}),
           'lee_input_sha256': digest(LEE), 'siraj_input_sha256': digest(OTHER),
           'output_sha256': digest(OUT), 'code_sha256': digest(Path(__file__))}
(HERE / 'cross_assay_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(receipt, indent=2))

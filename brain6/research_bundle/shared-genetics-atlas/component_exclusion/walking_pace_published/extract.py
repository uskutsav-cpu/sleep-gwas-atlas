"""Read selected published genetic correlations from SD8; do not edit workbook."""
import csv
import hashlib
import json
from pathlib import Path
import openpyxl

source = Path('/Volumes/Brain6Work/walking_pace_published_supplement.xlsx')
out = Path(__file__).resolve().parent
payload = source.read_bytes()
assert hashlib.md5(payload).hexdigest() == '3004d1fd880c636ed3fefd039dd4f6a7'
w = openpyxl.load_workbook(source, read_only=True, data_only=True)
s = w['SD8']
targets = {'Insomnia', 'Attention deficit hyperactivity disorder', 'Hand grip strength'}
rows = []
for cells in s.iter_rows(min_row=7, max_col=12):
    v = [c.value for c in cells]
    if v[1] in targets:
        rows.append(dict(sheet='SD8', row=cells[0].row, trait=v[1], study_n=v[2],
                         cohort=v[3], reference=v[4], pmid=v[5],
                         rg=v[6], se=v[7], p=v[8], bmi_adjusted_rg=v[9],
                         bmi_adjusted_se=v[10], bmi_adjusted_p=v[11]))
assert len(rows) == len(targets) and {x['trait'] for x in rows} == targets
assert all(x['p'] < 0.05 / 53 for x in rows)
assert all(x['bmi_adjusted_p'] < 0.05 / 53 for x in rows)
with (out / 'published_correlations.tsv').open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
    writer.writeheader(); writer.writerows(rows)
receipt = {'source_url': 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs42003-020-01357-7/MediaObjects/42003_2020_1357_MOESM4_ESM.xlsx',
           'source_file': str(source), 'source_bytes': len(payload),
           'source_md5_matches_europe_pmc_xml': hashlib.md5(payload).hexdigest(),
           'source_sha256': hashlib.sha256(payload).hexdigest(),
           'sheet': 'SD8', 'rows': [r['row'] for r in rows],
           'published_threshold_53_traits': 0.05 / 53,
           'output_sha256': hashlib.sha256((out / 'published_correlations.tsv').read_bytes()).hexdigest(),
           'note': 'Publication results; not a new LDSC analysis.'}
(out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(receipt, indent=2))

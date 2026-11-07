import csv,hashlib,json
from pathlib import Path
p=Path(__file__).resolve().parent; root=p.parents[1];receipt=json.loads((p/'receipt.json').read_text());coverage=json.loads((p.parent/'coverage_receipt.json').read_text());source=json.loads((root/'receipt.json').read_text())['sources'];h=receipt['source_sha256'];assert h['GSE244011_RAW.tar']==coverage['source_archive_sha256'] and h['unique_bcdat.tsv.gz']==coverage['barcode_map_sha256'];assert all(h[k]==source[k] for k in ('dnaMatrix.csv.gz','rnaMatrix.csv.gz'))
with (p/'counts.tsv').open(newline='') as f:r=list(csv.DictReader(f,delimiter='\t'))
assert len(r)==44
for s in receipt['summary']:
 q=[x for x in r if x['kind']==s['kind']];a=s['allele'];assert len(q)==22 and sum(int(x[a+'_raw']) for x in q)==s['raw_sum'] and sum(int(x[a+'_matrix']) for x in q)==s['matrix_sum'];assert abs(s['raw_to_matrix_sum_ratio']-s['raw_sum']/s['matrix_sum'])<1e-12
(p/'verification.json').write_text(json.dumps({'source_hashes_match_prior_receipts':True,'counts_rows':len(r),'summary_arithmetic_verified':True,'scope':'Not a validated sample crosswalk or full barcode reconstruction'},indent=2)+'\n');print('four source hashes and summarized counts verified')

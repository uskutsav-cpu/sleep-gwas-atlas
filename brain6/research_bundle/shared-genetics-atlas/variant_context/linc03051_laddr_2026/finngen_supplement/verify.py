#!/usr/bin/env python3
import csv,hashlib,json
from pathlib import Path
p=Path(__file__).resolve().parent;receipt=json.loads((p/'receipt.json').read_text())
def rows(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
g=rows(p/'gene_rows.tsv');m=rows(p/'matched_sleep_models.tsv')
assert len(g)==12 and all(x['gene_id']=='ENSG00000239268' and x['finngen_trait']=='F5_DEPRESSION_DYSTHYMIA' for x in g)
assert len(m)==2 and {x['tissue'] for x in m}=={'BRNHPP','BRNACC'}
assert all(float(x['finngen_clpp'])>.01 and x['feature']=='ENSG00000239268__PC2' for x in m)
for file,key in [('gene_rows.tsv','gene_rows_sha256'),('matched_sleep_models.tsv','matched_sha256')]:
 assert hashlib.sha256((p/file).read_bytes()).hexdigest()==receipt[key]
assert receipt['rows_read']==63741 and receipt['unique_reported_traits']==692
out={'matched_rows':len(m),'all_gene_rows':len(g),'rows_checked':receipt['rows_read'],
 'source_receipt_sha256':hashlib.sha256((p/'receipt.json').read_bytes()).hexdigest(),
 'note':'Verifies extracted data; not independent validation of publication or full testing coverage.'}
(p/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

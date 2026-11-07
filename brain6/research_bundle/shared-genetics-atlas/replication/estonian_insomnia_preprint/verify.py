#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
import openpyxl
p=Path(__file__).resolve().parent;r=json.loads((p/'receipt.json').read_text())
with (p/'matched_candidates.tsv').open() as f:matched=list(csv.DictReader(f,delimiter='\t'))
assert len(matched)==2 and {x['rsid'] for x in matched}=={'rs1093463','rs1989903'}
w=openpyxl.load_workbook(p/'media-2.xlsx',read_only=True,data_only=True)
rows=[x for x in w['ST5'].iter_rows(values_only=True) if x and x[0]=='Insomnia'];assert len(rows)==33
for x in matched:
 hit=[a for a in rows if a[3]==x['rsid'] and int(a[1])==int(x['chromosome']) and int(a[2])==int(x['position38']) and a[4]==x['allele1'] and a[5]==x['allele0']]
 assert len(hit)==1 and math.isclose(float(hit[0][14]),float(x['p']),rel_tol=1e-12) and x['is_one_of_six_leads']=='False'
assert hashlib.sha256((p/'matched_candidates.tsv').read_bytes()).hexdigest()==r['rows_sha256']
o={'matched_records_reproduced':2,'source_insomnia_rows':33,'matched_source_p_values':True,'none_six_leads':True,'limitation':'Checks selected table rows and matching, not source phenotype validity, cohort independence, genetic replication or statistical power.'}
(p/'verification.json').write_text(json.dumps(o,indent=2)+'\n');print(json.dumps(o,indent=2))

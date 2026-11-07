#!/usr/bin/env python3
import csv,hashlib,json,math,zipfile
from pathlib import Path
import openpyxl
p=Path(__file__).resolve().parent
with (p/'candidate_comparison.tsv').open() as f:r=list(csv.DictReader(f,delimiter='\t'))
assert len(r)==42 and len({x['rsid'] for x in r})==42
q=[x for x in r if x['allelic_status']=='unique_match'];assert len(q)==27
assert sum(float(x['allelic_fdr'])<.05 for x in q)==7
assert sum((float(x['allelic_fdr'])<.05)!=(float(x['v124_fdr'])<.05) for x in q)==0
assert sum((float(x['allelic_fdr'])<.05)!=(float(x['v134_fdr'])<.05) for x in q)==5
em=[x for x in q if x['allelic_class']=='emVar'];assert {x['rsid'] for x in em}=={'rs325485','rs6874138'}
assert all(x['element_status']=='unique_match' and x['element_class']=='MPRA-active' and float(x['allelic_fdr'])<.05 and float(x['element_fdr'])<.05 for x in em)
workbook=openpyxl.load_workbook(p/'table_s1.xlsx',read_only=True,data_only=True)
targets={x['rsid']:x for x in em};seen=set()
for row in workbook['Table S1B'].iter_rows(min_row=2,values_only=True):
 if row[3] in targets:
  x=targets[row[3]]; assert row[16]=='emVar' and math.isclose(row[15],float(x['allelic_fdr']),rel_tol=1e-12);seen.add(row[3])
assert seen==set(targets)
ld=json.loads((p/'ld_receipt.json').read_text())
assert len(ld['pairwise'])==4 and all(x['n_individuals']==503 and math.isclose(x['r']**2,x['r2'],abs_tol=1e-12) for x in ld['pairwise'])
receipt={'published_matched':27,'published_fdr_lt_005':7,'emvar_rsids':sorted(targets),'threshold_disagreements_v124':0,'threshold_disagreements_v134':5,'ld_pairs_arithmetic_checked':4,'sha256':{x:hashlib.sha256((p/x).read_bytes()).hexdigest() for x in ('table_s1.xlsx','candidate_comparison.tsv','ld_receipt.json')},'limitations':'Checks subset matching and arithmetic, not publisher QC, biological significance, independence, or underlying reference genotypes.'}
(p/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

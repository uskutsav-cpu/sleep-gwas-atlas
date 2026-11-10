#!/usr/bin/env python3
"""Extract author-reported SleepChart LDSC table, preserving source values."""
import csv, hashlib, json
from pathlib import Path
import openpyxl
src=Path('work/sleepchart_2026_supplement.xlsx')
out=Path('outputs/sleepchart_2026_author_ldsc.csv')
s=openpyxl.load_workbook(src,read_only=True,data_only=True)['SF5']
rows=list(s.values)
assert rows[2][:7]==('Duration type','DE','Dataset','gc_mean','gc_std','Z','P')
results=[list(r[:7]) for r in rows[3:] if r[0] in ('short_vs_normal','long_vs_normal') and r[1] is not None]
assert len(results)==1054 and len({(r[0],r[1],r[2]) for r in results})==1054
with out.open('w',newline='') as f:
 w=csv.writer(f); w.writerow(['duration_type','disease_endpoint','dataset','genetic_correlation','standard_error','z','p']); w.writerows(results)
receipt={'source':'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-026-10524-5/MediaObjects/41586_2026_10524_MOESM3_ESM.xlsx','sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'sheet':'SF5','columns':'A:G','rows':len(results),'output':str(out),'analysis':'Author-reported estimates, not new LDSC analysis'}
Path('outputs/sleepchart_2026_ldsc_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))

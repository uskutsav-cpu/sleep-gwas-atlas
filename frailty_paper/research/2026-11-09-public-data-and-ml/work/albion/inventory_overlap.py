#!/usr/bin/env python3
"""Aggregate publicly provided ALBION assay flags, without exporting identifiers."""
import hashlib,json,collections
from pathlib import Path
import openpyxl
p=Path(__file__).resolve().parent/'assay_inventory.xlsx'
w=openpyxl.load_workbook(p,read_only=True,data_only=True); it=w['Biomarkers'].iter_rows(values_only=True);header=next(it); rows=[r for r in it if r[0] is not None]
fields=['PAT administration','7 day Actigraphy','Imputed GWAS','Whole Genome Sequencing']; idx={name:header.index(name) for name in fields}
ids=collections.Counter(str(r[0]) for r in rows)
flag=lambda r,k:r[idx[k]]==1
stats={'inventory_rows_with_identifier':len(rows),'distinct_identifiers':len(ids),'identifiers_appearing_multiple_times':sum(c>1 for c in ids.values()),'flag_1_rows':{k:sum(flag(r,k) for r in rows) for k in fields},'rowwise_flag_1_overlap':{a+' & '+b:sum(flag(r,a) and flag(r,b) for r in rows) for a in fields[2:] for b in fields[:2]}}
result={'source':'https://albionstudy.github.io/ALBION_Biomarkers_JUL2026.xlsx','sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'sheet':'Biomarkers','flag_rule':'numeric 1 in both named columns in same row; missing and 0 not counted','statistics':stats,'limitations':'One identifier is duplicated; rows are assay availability flags rather than verified unique participants or usable observations. Unknown linkage, quality, access, consent, temporal alignment and score coverage.'}
out=Path(__file__).resolve().parents[2]/'outputs/stage6_albion_assay_overlap_receipt.json';out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

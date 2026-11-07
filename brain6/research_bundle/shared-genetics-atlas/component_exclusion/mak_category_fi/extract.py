#!/usr/bin/env python3
"""Read-only extraction from the published workbook and the versioned Figshare metadata."""
import csv,hashlib,json,pathlib
import openpyxl
p=pathlib.Path(__file__).resolve().parent;book=p/'source_supplement.xlsx';m=p/'figshare_v1_metadata.json'
w=openpyxl.load_workbook(book,read_only=True,data_only=True);s=w['ST4'];rows=[]
for row in range(12,23):
 label=s.cell(row,2).value;indicator=s.cell(row,3).value
 if indicator=='Yes':rows.append({'source_cell':f'B{row}:C{row}','item':label,'ukb_fi':'Yes'})
assert len(rows)==8 and {r['item'] for r in rows}.issuperset({'Fatigue','Insomnia'})
t=w['ST16'];title=t['A1'].value;header=t['Q2'].value
assert '8 mental wellbeing-related items' in header and 'A GWAS' in title
meta=json.loads(m.read_text());files=[f['name'] for f in meta['files']];assert len(files)==19
with (p/'excluded_items.tsv').open('w',newline='') as h:
 writer=csv.DictWriter(h,fieldnames=['source_cell','item','ukb_fi'],delimiter='\t');writer.writeheader();writer.writerows(rows)
receipt={'workbook_sha256':hashlib.sha256(book.read_bytes()).hexdigest(),'source_workbook_initial_path':'/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/mak2023_supplement.xlsx','supplement_sheet_header':header,'supplement_table_note':'ST16 A1 describes GWAS of 11 category-stripped FIs; values in ST16 are biomarker MR results, not per-SNP GWAS.','removed_item_count':len(rows),'nominal_residual_item_count':49-len(rows),'figshare_metadata_url':'https://api.figshare.com/v2/articles/24591573/versions/1','figshare_metadata_sha256':hashlib.sha256(m.read_bytes()).hexdigest(),'figshare_file_names':files,'caution':'Filename inventory does not prove absence in other repositories or that samples match; no genetic correlation was estimated.'}
(p/'extraction_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'n_items':len(rows),'filenames':len(files),'sha256':receipt['workbook_sha256']},indent=2))

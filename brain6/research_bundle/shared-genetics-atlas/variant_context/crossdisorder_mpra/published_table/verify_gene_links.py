#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
import openpyxl
p=Path(__file__).resolve().parent
r=json.loads((p/'gene_link_receipt.json').read_text())
w=openpyxl.load_workbook(p/'table_s4.xlsx',read_only=True,data_only=True)
targets={'rs325485','rs6874138'};found=[]
for sheet in ('Table S4A','Table S4B','Table S4C'):
 rows=w[sheet].iter_rows(values_only=True);header=next(rows);ix={v:i for i,v in enumerate(header)}
 for row in rows:
  if row[ix['RSID']] in targets:found.append((sheet,row[ix['RSID']],row[ix['CHR']],row[ix['BP']],row[ix['ENSG']] if 'ENSG' in ix else None,row[ix['Method']] if 'Method' in ix else None))
assert found==[('Table S4C','rs6874138',5,103899596,'ENSG00000112874','ABC')]
assert len(r['matched_rows'])==1 and not r['identifier_matches_with_coordinate_or_allele_mismatch']
receipt={'independently_matched_rows':len(found),'matching_gene_identifier':found[0][4],'method':found[0][5],'table_sha256':hashlib.sha256((p/'table_s4.xlsx').read_bytes()).hexdigest(),'receipt_sha256':hashlib.sha256((p/'gene_link_receipt.json').read_bytes()).hexdigest(),'limitation':'Checks table retrieval and row identity only, not correctness of gene mapping, gene nomenclature, association causality or model calibration.'}
(p/'gene_link_verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

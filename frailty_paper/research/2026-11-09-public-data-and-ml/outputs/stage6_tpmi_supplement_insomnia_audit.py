#!/usr/bin/env python3
"""Extract published insomnia/movement rows for a metadata and AUC audit."""
import csv, hashlib, json
from pathlib import Path
import openpyxl

root=Path(__file__).resolve().parent.parent
source=root/'work/stage6_tpmi_supplement.xlsx'
book=openpyxl.load_workbook(source,read_only=True,data_only=True)
selected=[]
for sheet in ('Table.S1','Table.S9','Table.S10','Table.S11','Table.S12','Table.S13'):
    for row in book[sheet].iter_rows(values_only=True):
        if str(row[0]) not in ('327.4','327.41','327.7'): continue
        if sheet=='Table.S1':
            selected.append(dict(sheet=sheet,phecode=row[0],phenotype=row[1],gwas_cases=row[12],gwas_controls=row[13],training_cases=row[14],training_controls=row[15],validation_cases=row[16],validation_controls=row[17]))
        else:
            offset=1 if sheet=='Table.S11' else 0
            selected.append(dict(sheet=sheet,phecode=row[0],phenotype=row[1],score_auc=row[2+offset],score_auc_se=row[3+offset],score_covariates_auc=row[6+offset],score_covariates_auc_se=row[7+offset]))
out=root/'outputs/stage6_tpmi_supplement_insomnia_audit.json'
out.write_text(json.dumps({'source':'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-025-09350-y/MediaObjects/41586_2025_9350_MOESM3_ESM.xlsx','sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'rows':selected},indent=2)+'\n')
print(f'Extracted {len(selected)} rows from {source.name}')

import openpyxl,csv,hashlib,json
from pathlib import Path
p=Path('work/stage6_gong2026_msk_supplement.xlsx');o=Path('outputs');w=openpyxl.load_workbook(p,read_only=True,data_only=True);out={}
for no,field,wanted in [(6,1,'rs12536335'),(7,2,'rs12536335'),(9,1,153),(14,2,'RBM6')]:
 s=w[f'Supplementary Data {no}'];rows=list(s.iter_rows(min_row=4,values_only=True));h=list(next(s.iter_rows(min_row=3,max_row=3,values_only=True)));matches=[r for r in rows if r[field]==wanted and (no!=9 or r[0]=='MSKdg')];name=f'stage6_gong2026_supplement_table{no}_selected.csv'
 with (o/name).open('w',newline='') as f:
  c=csv.writer(f);c.writerow(h);c.writerows(matches)
 out[str(no)]={'scanned_rows':len(rows),'matching_rows':len(matches),'file':name}
receipt={'supplement_url':'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-026-72164-7/MediaObjects/41467_2026_72164_MOESM2_ESM.xlsx','sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'selection':out,'warning':'Supplementary Data 14 contains significant rows only; no-row is not a negative result; factor samples overlap with project discovery.'};(o/'stage6_gong2026_source_table_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

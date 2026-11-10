import openpyxl,csv,hashlib,pathlib
src=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/foote2025_tables.xlsx'); w=openpyxl.load_workbook(src,read_only=True,data_only=True);s=w['Supplementary Table 44']; rows=[]
for i,row in enumerate(s.iter_rows(values_only=True),1):
 if row[0] in ('Multi-PRS','FI-PRS','FP-PRS') and row[1]=='Frailty Index':
  rows.append({'source_row':i,'score':row[0],'outcome':row[1],'cohort':row[2],'beta':row[3],'se':row[4],'ci_lower':row[5],'ci_upper':row[6],'p':row[7],'q':row[8],'r2_as_reported':row[9]})
assert len(rows)==9
p=pathlib.Path('outputs/foote_same_outcome_multiscore_comparison.csv')
with p.open('w',newline='') as f:
 x=csv.DictWriter(f,fieldnames=rows[0]);x.writeheader();x.writerows(rows)
print('workbook SHA256:',hashlib.sha256(src.read_bytes()).hexdigest());print('csv SHA256:',hashlib.sha256(p.read_bytes()).hexdigest());print(rows)

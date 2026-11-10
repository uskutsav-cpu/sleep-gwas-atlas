#!/usr/bin/env python3
import csv,hashlib,io,subprocess
from pathlib import Path
import openpyxl
root=Path(__file__).resolve().parent.parent;p=root/'work/frailpol_participants.xlsx'
if hashlib.md5(p.read_bytes()).hexdigest()!='c438dc375281f2fb1e0afe08be89aed3':raise ValueError('checksum')
s=openpyxl.load_workbook(p,read_only=True,data_only=True).active
rows=[r for r in list(s.values)[1:] if r[0] is not None and str(r[0]).isdigit()];lookup={str(r[0]):r for r in rows}
with (root/'work/frailpol_gait_three_class.csv').open(newline='') as f:gait=list(csv.DictReader(f))
if len(lookup)!=len(rows) or len(gait)!=668:raise ValueError('duplicate ID or unexpected count')
out=io.StringIO();writer=csv.writer(out);writer.writerow(['participant_id','eligible'])
for item in gait:
 r=lookup[item['participant_id']];valid=all(r[i] in (0,1) for i in range(9,14)) and isinstance(r[14],(int,float));writer.writerow([item['participant_id'],'TRUE' if valid and sum(r[9:14])==r[14] else 'FALSE'])
subprocess.run(['Rscript','work/frailpol_label_sensitivity.R'],cwd=root,input=out.getvalue(),text=True,check=True)

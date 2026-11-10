import openpyxl,math,csv,hashlib,pathlib
s=pathlib.Path('work/hfrs_2025_tables.xlsx');w=openpyxl.load_workbook(s,read_only=True,data_only=True);rows=[]
for i,row in enumerate(w['ST15'].iter_rows(values_only=True),1):
 if i<3 or row[2]!='Full sample':continue
 out={'table':'ST15','source_row':i,'outcome':row[0],'score':row[1],'population':row[2],'model':row[3],'beta':row[4],'OR_HR':row[5],'SE':row[6],'p_reported':row[7], 'minus_log10p_reported':row[8]}
 if row[3] in ('Logistic regression','Cox regression'):
  b,se=float(row[4]),float(row[6]);out['approx_wald_OR_HR_lower']=math.exp(b-1.95996398454*se);out['approx_wald_OR_HR_upper']=math.exp(b+1.95996398454*se);assert abs(math.exp(b)-float(row[5]))<1e-12
 rows.append(out)
p=pathlib.Path('outputs/hfrs_score_published_outcomes.csv'); keys=list(dict.fromkeys(k for row in rows for k in row))
with p.open('w',newline='') as f:wr=csv.DictWriter(f,fieldnames=keys);wr.writeheader();wr.writerows(rows)
print('source_sha256',hashlib.sha256(s.read_bytes()).hexdigest(),'output_sha256',hashlib.sha256(p.read_bytes()).hexdigest())
for x in rows:print(x)

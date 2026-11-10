import openpyxl,csv
p='/Volumes/Extreme SSD/Codex Archive/foote2025_tables.xlsx';w=openpyxl.load_workbook(p,read_only=True,data_only=True);out=[]
for number,cohort in [(52,'LBC1936'),(53,'ELSA'),(54,'PISA')]:
 for i,row in enumerate(w[f'Supplementary Table {number}'].iter_rows(values_only=True),1):
  if i<4:continue
  vals=[str(v) if v is not None else '' for v in row];joined=' '.join(vals).lower()
  if any(v in joined for v in ['sleep','insom','tired','fatigue','energy','exhaust']):
   out.append(dict(cohort=cohort,source_table=number,source_row=i,source_column_1=vals[0],source_column_2=vals[1],source_column_3=vals[2],caution='Keyword identified item; not necessarily an insomnia diagnosis or independently verified recoding'))
with open('outputs/foote_fi_evaluation_sleep_items.csv','w',newline='') as f:
 q=csv.DictWriter(f,fieldnames=list(out[0]));q.writeheader();q.writerows(out)
print('total',len(out),'ELSA',sum(x['cohort']=='ELSA' for x in out),'PISA',sum(x['cohort']=='PISA' for x in out),'LBC1936',sum(x['cohort']=='LBC1936' for x in out))

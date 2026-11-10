import openpyxl,csv
w=openpyxl.load_workbook('/Volumes/Extreme SSD/Codex Archive/foote2025_tables.xlsx',read_only=True,data_only=True);s=w['Supplementary Table 44'];out=[]
for no,row in enumerate(s.iter_rows(values_only=True),1):
 if row[0] in ('FI-PRS','FP-PRS'):
  out.append({'source_row':no,'predictor':row[0],'outcome':row[1],'cohort':row[2],'estimate':row[3],'se':row[4],'ci_lower':row[5],'ci_upper':row[6],'p_value':row[7],'fdr_q':row[8],'variance_explained_reported':row[9],'source':'Foote 2025 Supplementary Table 44','caution':'Reported values; predictor labels interpreted via paper; independent paired comparison not supplied'})
with open('outputs/foote_supplement_44_fi_ffs_score_results.csv','w',newline='') as f:
 q=csv.DictWriter(f,fieldnames=list(out[0]));q.writeheader();q.writerows(out)
print('extracted',len(out),'rows')

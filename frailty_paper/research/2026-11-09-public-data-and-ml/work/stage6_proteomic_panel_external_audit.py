import io,zipfile,json,openpyxl
from pathlib import Path
with zipfile.ZipFile('work/pfi_supplement.zip') as z:b=z.read('ACEL-24-e70144-s001.xlsx')
s=openpyxl.load_workbook(io.BytesIO(b),data_only=True,read_only=True)['Supplementary Table 8']
r=[row for row in list(s.values)[3:] if isinstance(row[5],(float,int))]
assert len(r)==25
p=.05/25
same=[x for x in r if x[5]*x[7]>0 and x[5]*x[9]>0]
both=[x for x in same if x[8]<p and x[10]<p]
o={'n':25,'same_direction_across_3':len(same),'both_external_bonferroni_0_05_over_25':len(both),'exceptions_direction':[x[2] for x in r if x not in same],'passing_both_external':[x[2] for x in both],'threshold':p,'note':'Post-hoc descriptive screen in a preselected 25-protein panel; measured-protein associations, not genetic or sleep results. External results differ in power and covariates.'}
Path('outputs/stage6_proteomic_external_panel_audit.json').write_text(json.dumps(o,indent=2)+'\n');print(json.dumps(o,indent=2))

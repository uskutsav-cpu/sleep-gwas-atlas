import hashlib,io,json,subprocess
from pathlib import Path
import openpyxl
p=Path(__file__).resolve().parent
u='https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-59383-0/MediaObjects/41467_2025_59383_MOESM3_ESM.xlsx'
b=subprocess.check_output(['curl','-fsSL','--max-time','50',u]);h=hashlib.sha256(b).hexdigest();assert h=='3bb710a7c1871a79936c0ac35ba283865344aca667c4b770718d46d14a8fe8b1'
w=openpyxl.load_workbook(io.BytesIO(b),read_only=True,data_only=True);s=w['Supplementary Data 2'];found=[]
for row in s.iter_rows(values_only=True):
 if len(row)>7 and str(row[7]).strip().lower()=='physical':found.append({'phenotype':row[1],'minimum_years':row[2],'maximum_years':row[3],'mean_years':row[4],'complete_followup_n':row[6]})
assert {x['phenotype'] for x in found}=={'Forced expiratory volume','Hand grip strength','Height'}
(p/'supplement_receipt.json').write_text(json.dumps({'url':u,'bytes':len(b),'sha256':h,'sheet':'Supplementary Data 2','physical_rows':found,'scope':'The paper states baseline composite uses the same selected aging phenotypes; table alone does not verify source genotype rows or individual score calculations.'},indent=2)+'\n');print(json.dumps(found,indent=2))

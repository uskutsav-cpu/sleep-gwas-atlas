import json,csv,hashlib
from pathlib import Path
selected={'ICD_codes':{'G47','F51','R53'},'Phecodes':{'327','798'}};out=[];receipts={}
for score in ['Fried_Frailty','Insomnia']:
 for typ,ids in selected.items():
  p=Path(f'work/polygenie_{score}_{typ}.json');j=json.loads(p.read_text())['response']['prs-table-content']['data'];receipts[f'{score}/{typ}']={'rows':len(j),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
  for row in j:
   if str(row['Code']) in ids:out.append({'score_label':score,'outcome_type':typ,**row})
assert len(out)==10
with open('outputs/stage6_polygenie_sleep_frailty_results.csv','w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
Path('outputs/stage6_polygenie_sleep_frailty_receipt.json').write_text(json.dumps({'source':'https://polygenie.igtp.cat/','comparison':'top PRS decile versus bottom decile','files':receipts,'limitations':'post-hoc broad scores and diagnoses; controlled individual genotypes; cohort and relative overlap unverified'},indent=2)+'\n')
print('exported',len(out),'rows')

import csv,hashlib,json
from pathlib import Path
import openpyxl
p=Path('work/koges/nam2022_table_s1.xlsx')
s=openpyxl.load_workbook(p,read_only=True,data_only=True).active
spans=[('disease',6,19),('biomarker',24,54),('diet',59,81),('alcohol',86,93)]
rows=[]
for category,a,b in spans:
 for row in s.iter_rows(min_row=a,max_row=b,values_only=True):
  assert row[0] and row[1]
  rows.append({'category':category,'phenotype':row[0],'type':row[1]})
assert len(rows)==76
out=Path('outputs/stage6_koges_nam2022_76_phenotypes.csv')
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['category','phenotype','type']);w.writeheader();w.writerows(rows)
terms=('insomnia','sleep','frail','grip','handgrip')
hits=[r for r in rows if any(t in r['phenotype'].lower() for t in terms)]
receipt={'source_url':'https://ars.els-cdn.com/content/image/1-s2.0-S2666979X22001380-mmc2.xlsx','source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'table':'S1. Basic information','rows':len(rows),'category_counts':{x:sum(r['category']==x for r in rows) for x,_,_ in spans},'search_terms':terms,'matching_phenotype_labels':hits,'note':'Search only declared phenotype labels and not a study-wide census or individual level validation.'}
Path('outputs/stage6_koges_nam2022_76_phenotypes_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))

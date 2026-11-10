import csv,json,hashlib
from pathlib import Path
import openpyxl
p=Path('work/kunorozva2025/supplement.xlsx')
files={'score15':'outputs/unvalidated_15_variant_score_prototype.csv','count22':'outputs/unvalidated_22_shared_direction_score.csv'}
sets={k:{r['variant'] for r in csv.DictReader(open(v))} for k,v in files.items()};targets=set.union(*sets.values())
assert len(sets['score15'])==15 and len(sets['count22'])==22 and len(targets)==33
w=openpyxl.load_workbook(p,read_only=True,data_only=True);out=[];counts={}
for sheet,offset,index in [('ST3',3,2),('ST5',2,3)]:
 rows=w[sheet].iter_rows(values_only=True);header=None;seen=0
 for i,row in enumerate(rows,1):
  if i==offset:header=row;continue
  if i<=offset or len(row)<=index or not row[index]:continue
  seen+=1;vid=str(row[index])
  if vid not in targets:continue
  cols=15
  out.append({'sheet':sheet,'source_row':i,'variant':vid,'in_score15':vid in sets['score15'],'in_count22':vid in sets['count22'],'source_fields_json':json.dumps(dict(zip([str(x) for x in header[:cols]],row[:cols])),default=str,ensure_ascii=False)})
 counts[sheet]=seen
assert sum(r['sheet']=='ST3' for r in out)==5 and sum(r['sheet']=='ST5' for r in out)==3
output=Path('outputs/stage6_kunorozva2025_candidate_selected_rows.csv')
with output.open('w',newline='') as f:
 x=csv.DictWriter(f,fieldnames=out[0].keys());x.writeheader();x.writerows(out)
receipt={'source_url':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12668095/supplementaryFiles','source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'n_unique_score_variants':33,'n_rows_by_sheet':counts,'n_matched_rows':len(out),'matched_rsids':sorted({r['variant'] for r in out}),'search':'Exact rsID in selected meta-analysis lead and Estonian follow-up tables only','warning':'No complete GWAS, allele semantics unverified, phenotype and cohort differences; no score validation.'}
Path('outputs/stage6_kunorozva2025_candidate_selected_rows_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

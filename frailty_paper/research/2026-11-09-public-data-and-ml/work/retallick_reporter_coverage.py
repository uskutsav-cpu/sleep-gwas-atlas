import csv,openpyxl,re,json,hashlib
from pathlib import Path
p=Path('outputs');book=Path('work/retallick_neural/published_supp2.xlsx');ids={r['variant'] for f in ['integrated_15_candidate_review.csv','expanded_22_candidate_risk_gene_table.csv'] for r in csv.DictReader((p/f).open())};assert len(ids)==33
w=openpyxl.load_workbook(book,read_only=True,data_only=True);out={}
for sn in ['SI Data 2.6','SI Data 2.8','SI Data 2.9','SI_Data 2.12','SI_Data 2.13','SI_Data 2.14','SI_Data 2.15','SI_Data 2.16']:
 matches={}
 for i,row in enumerate(w[sn].iter_rows(values_only=True),1):
  inrow=set()
  for value in row:
   if isinstance(value,str):inrow.update(set(re.findall(r'(?<![A-Za-z0-9])rs\d+(?![A-Za-z0-9])',value))&ids)
  for v in inrow:matches.setdefault(v,[]).append(i)
 out[sn]={'n_matched_variants':len(matches),'matched_sheet_rows':matches}
r={'publication':'https://doi.org/10.1038/s41467-026-77114-x','supplement_url':'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-026-77114-x/MediaObjects/41467_2026_77114_MOESM5_ESM.xlsx','published_workbook_sha256':hashlib.sha256(book.read_bytes()).hexdigest(),'project_candidate_count':len(ids),'matching':'Exact rsID text in each worksheet; no build/allele proof and no inference from absence.','worksheets':out}
(p/'stage6_retallick_candidate_coverage.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

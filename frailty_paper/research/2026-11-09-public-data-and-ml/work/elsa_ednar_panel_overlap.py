import csv,hashlib,json,re
from pathlib import Path
import openpyxl
p=Path('work/elsa_ednar_frailty_supplement.xlsx'); ws=openpyxl.load_workbook(p,read_only=True,data_only=True)['Sheet1']
header=next(ws.iter_rows(min_row=13,max_row=13,values_only=True))
if list(header[:6])!=['CHR','SNP','BP','Allele','BETA','P']:raise ValueError('unexpected sheet layout')
rows=[r for r in ws.iter_rows(min_row=14,values_only=True) if isinstance(r[1],str) and re.fullmatch(r'rs\d+',r[1])]
ids=[r[1] for r in rows]
if len(ids)!=590 or len(set(ids))!=590:raise ValueError('unexpected table size or duplicate identifiers')
result={'article_doi':'10.1007/s40520-015-0419-z','supplement':'https://static-content.springer.com/esm/art%3A10.1007%2Fs40520-015-0419-z/MediaObjects/40520_2015_419_MOESM2_ESM.xlsx','supplement_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'table_rows':len(ids),'counts':{},'limits':'Exact rsID intersection with 590 analyzed records only; not the complete 620-marker panel, not an imputed dataset, no coordinate, allele, genotype, LD-proxy or population/relative-overlap confirmation. No patient-level score test.'}
for name in ['unvalidated_15_variant_score_prototype.csv','unvalidated_22_shared_direction_score.csv']:
 path=Path('outputs')/name
 with path.open(newline='') as f: score={r['variant'] for r in csv.DictReader(f)}
 result['counts'][name]={'score_n':len(score),'exact_rsid_overlap':len(score.intersection(ids)),'score_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
Path('outputs/stage6_elsa_590_marker_overlap.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

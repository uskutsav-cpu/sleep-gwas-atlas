#!/usr/bin/env python3
"""Compare reported lead identifiers with the two project score lists."""
import csv,hashlib,json
from pathlib import Path
import openpyxl
src=Path('work/sleepchart_2026_figures.xlsx')
s=openpyxl.load_workbook(src,read_only=True,data_only=True)['SD2']
it=s.iter_rows(values_only=True)
next(it)
h=next(it)
assert h[:6]==('Top Lead SNP','Chromosome','Position','P-value','Phenotype','Cytogenetic region')
rows=[]
for r in it:
 if r[0] is None: break
 rows.append(dict(zip(('id','chromosome','position_source_build_not_confirmed','p','phenotype','cytogenetic_region'),r[:6])))
assert len(rows)==8
lists={}
for label,path in [('15','outputs/unvalidated_15_variant_score_prototype.csv'),('22','outputs/unvalidated_22_shared_direction_score.csv')]:
 with open(path,newline='') as f: lists[label]={r['variant'] for r in csv.DictReader(f)}
for r in rows:
 r['exact_score15_rsid_match']=r['id'] in lists['15']; r['exact_score22_rsid_match']=r['id'] in lists['22']
result={'source':'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-026-10524-5/MediaObjects/41586_2026_10524_MOESM5_ESM.xlsx','source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'source_sheet':'SD2 first 8 lead variants','lead_count':len(rows),'exact_matches_15':sum(r['exact_score15_rsid_match'] for r in rows),'exact_matches_22':sum(r['exact_score22_rsid_match'] for r in rows),'warning':'Selected lead variants only; rsID identity comparison; nonmatch is not evidence of absent association or locus overlap.', 'rows':rows}
Path('outputs/sleepchart_2026_lead_variant_check.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))

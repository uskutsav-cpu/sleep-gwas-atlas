#!/usr/bin/env python3
"""Descriptive literal rsID overlap only; no genome build or LD inference."""
import csv,hashlib,json
from pathlib import Path
import openpyxl
p=Path(__file__).resolve().parent;root=p.parent.parent
paths={'paper':p/'source_S3.xlsx','brain':root/'variant_context/brain6_fi_ffs_screen/results.tsv','placo':root/'placo_diagnostic/completed_audit/candidate_positions.tsv'}
s=openpyxl.load_workbook(paths['paper'],read_only=True,data_only=True).active
rows=list(s.values);h=list(rows[0]);paper=[dict(zip(h,r)) for r in rows[1:] if r[1] is not None]
def read(x):
 with x.open() as f:return list(csv.DictReader(f,delimiter='\t'))
brain={r['variant'] for r in read(paths['brain'])};placo={r['SNP'] for r in read(paths['placo'])};rs={r['omni_name'] for r in paper if str(r['omni_name']).startswith('rs')}
res={'paper_rows':len(paper),'paper_literal_rsids':len(rs),'paper_other_identifiers':len(paper)-len(rs),'brain_variant_count':len(brain),'placo_unique_variant_count':len(placo),'literal_overlap_brain':sorted(rs&brain),'literal_overlap_placo':sorted(rs&placo),'input_sha256':{k:hashlib.sha256(x.read_bytes()).hexdigest() for k,x in paths.items()},'limitations':'Literal rsID only; unconverted kgp IDs, historical aliasing, linkage disequilibrium and regional overlaps not assessed; sample selection differs.'}
(p/'literal_id_check.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))

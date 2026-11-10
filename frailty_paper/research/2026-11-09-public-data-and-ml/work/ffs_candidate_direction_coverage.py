import csv, hashlib, json
from pathlib import Path
p=Path('outputs')
files=['ffs_cojo_effect_allele_audit.csv','unvalidated_15_variant_score_prototype.csv','unvalidated_22_shared_direction_score.csv']
def rows(f):
 with (p/f).open(newline='') as h:return list(csv.DictReader(h))
lead=rows(files[0]); scores=[rows(x) for x in files[1:]]
assert len(lead)==37 and all(len({r['variant'] for r in a})==len(a) for a in [lead]+scores)
l={r['variant'] for r in lead}
r={'match_definition':'exact variant identifier, not coordinate or LD','inputs':{f:hashlib.sha256((p/f).read_bytes()).hexdigest() for f in files},'lead_count':len(l),'scores':{f:{'count':len(a),'matching_identifiers':sorted(l.intersection(r['variant'] for r in a))} for f,a in zip(files[1:],scores)}}
(p/'stage6_ffs_candidate_direction_coverage.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

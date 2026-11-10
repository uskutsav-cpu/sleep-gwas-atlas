import json,math
from pathlib import Path
rows=json.loads(Path('work/mdd_symptoms/candidate_hits_private.json').read_text()); assert len(rows)==33
out={}
for trait in ['insomnia','fatigue']:
 for r in rows:assert len(r['hits'][trait])==1
 p=[float(r['hits'][trait][0]['p']) for r in rows]
 out[trait]={'exact_coordinate_and_pair_matches':len(p),'p_lt_0_05':sum(x<.05 for x in p),'p_lt_0_05_div_66':sum(x<.05/66 for x in p),'minimum_nominal_p':min(p)}
Path('outputs/stage6_mdd_clinical_candidate_summary.json').write_text(json.dumps({'source':'https://doi.org/10.6084/m9.figshare.22745573','traits':out,'threshold_note':'post hoc descriptive family of 33 markers × two traits; does not correct all prior screens','checksums_source':'Figshare API version 4 compressed-file computed_md5','limitations':'case-enriched symptom contrasts, not diagnosed primary insomnia, FI, FFS, or score prediction; cohorts and relatives overlap unverified'},indent=2)+'\n');print(json.dumps(out,indent=2))

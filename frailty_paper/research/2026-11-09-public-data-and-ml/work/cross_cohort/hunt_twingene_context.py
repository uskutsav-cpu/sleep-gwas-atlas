#!/usr/bin/env python3
"""Join published HUNT insomnia and TwinGene FI evidence for score-15 variants."""
import csv,hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[2]
paths={'hunt':root/'outputs/stage6_lane2019_hunt_partners_candidate_lookup.csv','twin':root/'outputs/candidate_twingene_fi_cohort_lookup.csv','score':root/'outputs/unvalidated_15_variant_score_prototype.csv'}
def load(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
h={r['variant']:r for r in load(paths['hunt']) if r['table_match']=='yes' and r['score_15_counted_other']};t={r['variant']:r for r in load(paths['twin'])};s={r['variant']:r for r in load(paths['score'])}; rows=[]
for v in sorted(h):
 a,b,c=h[v],t[v],s[v]
 assert a['chr_pos_source']==f"{c['chr']}:{c['position_GRCh37_inferred']}"
 assert b['matching_rows']=='1' and {a['source_effect_other'].split('/')[0],a['source_effect_other'].split('/')[1]}=={b['source_effect_allele'],b['source_other_allele']}
 rows.append({'variant':v,'hunt_effect_other':a['source_effect_other'],'hunt_or_ci':a['hunt_or_ci'],'hunt_insomnia_p':a['hunt_p'],'partners_insomnia_p':a['partners_p'],'twingene_fi_p':b['TG_p'],'ukb_fi_p':b['UKB_p'],'meta_fi_p':c['FI_p'],'ukb_ffs_p':c['FFS_p'],'notes':'Association p-values only; cohorts, ascertainment and selection differ; no harmonized joint effect or score validation.'})
out=root/'outputs/stage6_hunt_twingene_fi_candidate_context.csv'
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
receipt={'inputs':{k:{'path':str(p.relative_to(root)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for k,p in paths.items()},'output_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'rows':len(rows)}
(root/'outputs/stage6_hunt_twingene_fi_candidate_context_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))

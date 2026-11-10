import csv,json,pathlib
with open('outputs/stage6_zenodo_sarcopenia_candidates.csv') as f:
 rows=list(csv.DictReader(f))
assert len({r['rsid'] for r in rows})==len(rows)
lookup={r['rsid']:r for r in rows}; summary={}
for label,path in [('prototype15','outputs/unvalidated_15_variant_score_prototype.csv'),('prototype22','outputs/expanded_22_candidate_risk_gene_table.csv')]:
 with open(path) as f: items=list(csv.DictReader(f))
 matched=[];absent=[];mismatch=[]
 for row in items:
  name=row['variant'];src=lookup.get(name)
  if src is None:absent.append(name);continue
  c=row.get('chr',row.get('GRCh37_chr_assumed'));p=row.get('position_GRCh37_inferred',row.get('GRCh37_pos_assumed'))
  if str(c)==src['chromosome'] and str(p)==src['base_pair_location'] and {row['counted_allele'],row['other_allele']}=={src['effect_allele'],src['other_allele']}:matched.append(name)
  else:mismatch.append(name)
 summary[label]={'n':len(items),'id_position_pair_match':len(matched),'missing_rsid':absent,'mismatch':mismatch}
pathlib.Path('outputs/stage6_zenodo_sarcopenia_match_audit.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

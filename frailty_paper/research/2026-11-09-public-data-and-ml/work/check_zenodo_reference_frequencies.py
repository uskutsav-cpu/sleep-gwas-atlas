import csv,json,hashlib,pathlib
source='outputs/stage6_zenodo_sarcopenia_candidates.csv';reference='outputs/22_score_reference_group_allele_frequencies.csv';score='outputs/unvalidated_22_variant_score_prototype.csv'
src={r['rsid']:r for r in csv.DictReader(open(source))};ref={r['variant']:r for r in csv.DictReader(open(reference))}; weights={r['variant']:r for r in csv.DictReader(open(score))}
rows=[]
for v,w in weights.items():
 s=src[v];r=ref[v];assert {s['effect_allele'],s['other_allele']}=={w['counted_allele'],w['other_allele']} and r['counted_allele']==w['counted_allele']
 freq=float(r['reference_EUR_frequency']);ref_effect=freq if s['effect_allele']==w['counted_allele'] else 1-freq;reported=float(s['effect_allele_frequency']);delta=abs(reported-ref_effect)
 rows.append({'variant':v,'source_effect_allele':s['effect_allele'],'reported_eaf':reported,'reference_eur_effect_allele_frequency':ref_effect,'absolute_difference':delta,'within_half_rounding_unit_4dp':delta<=0.0000500000001})
with open('outputs/stage6_zenodo_eaf_reference_check.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
result={'n':len(rows),'within_0_00005':sum(r['within_half_rounding_unit_4dp'] for r in rows),'max_absolute_difference':max(r['absolute_difference'] for r in rows),'min_absolute_difference':min(r['absolute_difference'] for r in rows),'reference':'1000 Genomes phase 3, EUR 503 complete-call individuals, 1006 allele copies per variant','interpretation':'compatible with rounded reference frequencies, not proof of origin; 22 markers selected, source frequencies not established as study-derived','source_sha256':hashlib.sha256(open(source,'rb').read()).hexdigest(),'reference_sha256':hashlib.sha256(open(reference,'rb').read()).hexdigest()}
pathlib.Path('outputs/stage6_zenodo_eaf_reference_check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

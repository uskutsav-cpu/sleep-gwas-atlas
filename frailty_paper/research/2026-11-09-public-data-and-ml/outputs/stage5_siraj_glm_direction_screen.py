#!/usr/bin/env python3
"""Add descriptive source-reported ALT-vs-REF and model direction comparisons."""
import csv,collections
src='outputs/stage5_siraj2026_allcell_count_orientation.csv';dst='outputs/stage5_siraj2026_glm_direction_screen.csv'
with open(src) as f:rows=list(csv.DictReader(f));fields=list(rows[0])
extra=['glm_reported_direction','model_direction','descriptive_sign_comparison','reliability_scope']
with open(dst,'w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields+extra);w.writeheader()
 for x in rows:
  try:
   a=float(x['Log2Skew']);b=float(x['model_DNASE_raw']);ar='ALT higher' if a>0 else ('ALT lower' if a<0 else 'zero');br='ALT higher' if b>0 else ('ALT lower' if b<0 else 'zero');comp='same' if a*b>0 else ('opposite' if a*b<0 else 'zero')
  except (ValueError,TypeError):ar=br=comp='unavailable'
  x.update(glm_reported_direction=ar,model_direction=br,descriptive_sign_comparison=comp,reliability_scope='Study active-and-skew criterion met' if x['emvar_screen']=='TRUE' else 'Not classified as active-and-skew; sign is not significant evidence of reporter allelic effect');w.writerow(x)
print('groups',len(rows),'variants',len({r['variant_GRCh38'] for r in rows}))

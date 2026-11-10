#!/usr/bin/env python3
"""Cross-check reported skew against author-labeled aggregate ref/alt counts."""
import csv,math,collections,argparse
p=argparse.ArgumentParser();p.add_argument('--selected-counts',required=True);a=p.parse_args()
with open(a.selected_counts) as f:r=list(csv.DictReader(f))
with open('outputs/stage5_siraj2026_allcell_reporter_model_comparison.csv') as f:m=list(csv.DictReader(f))
model=collections.defaultdict(list)
for x in m:
 if x['source']=='GTEx' and x['cell_type']=='SKNSH':model[x['variant']].append(x)
groups={}
for x in r:
 if x['emvar_screen']!='TRUE':continue
 k=(x['variant'],x['cell_type'],x['library'],x['ID']);groups.setdefault(k,x)
fields=['variant_GRCh38','variant_GRCh37','library','mean_Plasmid_ref','mean_Plasmid_alt','mean_RNA_ref','mean_RNA_alt','aggregate_log2_alt_over_ref','reported_Log2Skew','count_based_vs_reported_sign','model_DNASE_alt_minus_ref','count_based_vs_model_sign','interpretation']
with open('outputs/stage5_siraj2026_five_variant_direction_check.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
 for k,x in groups.items():
  candidates=[y for y in model[x['variant']] if y['library']==x['library'] and y['emvar_screen']=='TRUE'];assert candidates
  g=candidates[0];pr,pa,rr,ra=[float(x[c]) for c in ['mean_Plasmid_ref','mean_Plasmid_alt','mean_RNA_ref','mean_RNA_alt']]
  assert min(pr,pa,rr,ra)>0
  delta=math.log2((ra/pa)/(rr/pr));reported=float(x['Log2Skew']);dnase=float(g['model_DNASE_raw'])
  w.writerow({'variant_GRCh38':g['variant_GRCh38'],'variant_GRCh37':x['variant'],'library':x['library'],**{c:x[c] for c in ['mean_Plasmid_ref','mean_Plasmid_alt','mean_RNA_ref','mean_RNA_alt']},'aggregate_log2_alt_over_ref':delta,'reported_Log2Skew':reported,'count_based_vs_reported_sign':'same' if delta*reported>0 else 'different_or_zero','model_DNASE_alt_minus_ref':dnase,'count_based_vs_model_sign':'same' if delta*dnase>0 else 'different_or_zero','interpretation':'Aggregate ratios are a descriptive cross-check, not the GLM effect or a paired uncertainty estimate; reporter RNA and DNase are distinct assays, and model training independence is unverified.'})
print('unique selected variant/library rows',len(groups))

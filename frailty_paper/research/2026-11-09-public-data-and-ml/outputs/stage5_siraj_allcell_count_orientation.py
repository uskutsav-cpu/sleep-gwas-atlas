#!/usr/bin/env python3
"""Cross-check allele-labeled reporter count means, not infer assay equivalence."""
import argparse,csv,collections,math
p=argparse.ArgumentParser();p.add_argument('--gtex',required=True);p.add_argument('--traits',required=True);a=p.parse_args()
fields=['mean_Plasmid_ref','mean_Plasmid_alt','mean_RNA_ref','mean_RNA_alt'];cols=['variant','cell_type','library','best_library','type','log2FC','logPadj_BF','Log2Skew','Skew_logPadj','active_screen','emvar_screen'];index=collections.defaultdict(set)
for path,source in [(a.gtex,'GTEx'),(a.traits,'complex_trait')]:
 with open(path) as f:
  for x in csv.DictReader(f):index[(source,*[x[k] for k in cols])].add(tuple(x.get(k,'') for k in fields))
with open('outputs/stage5_siraj2026_allcell_reporter_model_comparison.csv') as f:groups=list(csv.DictReader(f))
extra=fields+['aggregate_log2_alt_over_ref','count_vs_glm_direction','aggregate_vs_model_direction','orientation_status']
with open('outputs/stage5_siraj2026_allcell_count_orientation.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(groups[0])+extra);w.writeheader()
 for x in groups:
  counts=index.get((x['source'],*[x[k] for k in cols]),set());values=next(iter(counts)) if len(counts)==1 else ('',)*4
  row={**x,**dict(zip(fields,values)),'aggregate_log2_alt_over_ref':'','count_vs_glm_direction':'not assessable','aggregate_vs_model_direction':'not assessable','orientation_status':'ambiguous count records' if len(counts)>1 else ('missing counts' if not counts else 'positive labeled counts')}
  try:
   pr,pa,rr,ra=(float(v) for v in values);skew=float(x['Log2Skew']);dnase=float(x['model_DNASE_raw']);assert min(pr,pa,rr,ra)>0
   delta=math.log2((ra/pa)/(rr/pr));row['aggregate_log2_alt_over_ref']=delta
   if delta and skew:row['count_vs_glm_direction']='same' if delta*skew>0 else 'different'
   if delta and dnase:row['aggregate_vs_model_direction']='same' if delta*dnase>0 else 'different'
  except (ValueError,ZeroDivisionError,AssertionError):row['orientation_status']='missing or nonpositive counts'
  w.writerow(row)
print('groups',len(groups),'matched source keys',len(index))

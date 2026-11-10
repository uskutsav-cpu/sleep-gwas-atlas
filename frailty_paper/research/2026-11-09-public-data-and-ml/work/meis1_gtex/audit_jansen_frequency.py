import csv,gzip,json,math,hashlib
from pathlib import Path
p=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/data/raw/insomnia.txt.gz');rows=list(csv.DictReader(open('work/meis1_gtex/meis1_nfe_aligned.csv')))
idx={r['rs']:r for r in csv.DictReader(open('outputs/stage6_meis1_nfe_variant_index.csv'))}; need={(int(idx[r['id']]['pos37']),tuple(sorted([idx[r['id']]['ref37'],idx[r['id']]['alt37']]))):r for r in rows}; assert len(need)==len(rows); hits={}; n=0
with gzip.open(p,'rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  n+=1
  if x['CHR']!='2':continue
  try:k=(int(x['BP']),tuple(sorted([x['A1'],x['A2']])))
  except ValueError:continue
  if k in need:hits.setdefault(k,[]).append(x)
out=[]
for k,r in need.items():
 found=hits.get(k,[])
 if len(found)!=1:out.append({'id':r['id'],'status':'unmatched_or_duplicate'});continue
 x=found[0]
 try:maf=float(x['MAF']);nfe=float(r['nfe_alt_af']);refmaf=min(nfe,1-nfe);diff=abs(maf-refmaf)
 except ValueError:out.append({'id':r['id'],'status':'invalid_frequency'});continue
 out.append({'id':r['id'],'rsid':r['rsid'],'status':'matched','jansen_maf':maf,'nfe_maf':refmaf,'absolute_maf_difference':diff,'reported_unique_id':x['UNIQUE_ID']})
valid=[x for x in out if x['status']=='matched'];summary={'source':str(p),'rows_scanned':n,'variants':len(rows),'matched':len(valid),'unmatched':len(rows)-len(valid),'max_abs_maf_difference':max((x['absolute_maf_difference'] for x in valid),default=None),'n_difference_gt_0_2':sum(x['absolute_maf_difference']>0.2 for x in valid),'note':'Compared minor allele frequency only; MAF does not resolve effect-allele or strand orientation. External reference, source multi-allelic representations and sample overlap remain unverified.'}
Path('outputs/stage6_meis1_jansen_nfe_maf_audit.json').write_text(json.dumps({'summary':summary,'rows':out},indent=2)+'\n');print(json.dumps(summary,indent=2))

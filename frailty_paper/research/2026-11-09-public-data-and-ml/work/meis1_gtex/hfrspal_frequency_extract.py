import csv,json,hashlib
from pathlib import Path
p=Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw/GCST90707291.tsv');q={r['id']:r for r in csv.DictReader(open('work/meis1_gtex/meis1_hfrspal_input.csv'))};found={};n=0
with p.open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  n+=1
  if r['chromosome']!='2':continue
  for key in (f"2:{r['base_pair_location']}:{r['other_allele']}:{r['effect_allele']}",f"2:{r['base_pair_location']}:{r['effect_allele']}:{r['other_allele']}"):
   if key in q:found.setdefault(key,[]).append(r)
out=[]
for k,x in q.items():
 arr=found.get(k,[])
 if len(arr)!=1:continue
 r=arr[0]
 try:af=float(r['effect_allele_frequency'])
 except ValueError:continue
 if r['effect_allele']!=x['alt']:af=1-af
 out.append({'id':k,'hfrs_alt_af':af})
with open('work/meis1_gtex/meis1_hfrspal_frequencies.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=['id','hfrs_alt_af']);w.writeheader();w.writerows(out)
print(json.dumps({'rows_scanned':n,'variants':len(q),'unique_frequencies':len(out)}))

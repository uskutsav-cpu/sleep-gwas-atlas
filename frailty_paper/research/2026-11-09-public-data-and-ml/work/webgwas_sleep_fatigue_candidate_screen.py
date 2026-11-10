import csv, io, zipfile, hashlib, json, math
from pathlib import Path
archive=Path('work/webgwas_sleep_fatigue_20261029.zip')
files=['outputs/unvalidated_15_variant_score_prototype.csv','outputs/unvalidated_22_variant_score_prototype.csv']
ids={}
for fn in files:
 with open(fn,newline='') as f:
  for row in csv.DictReader(f): ids.setdefault(row['variant'],[]).append(Path(fn).name)
with zipfile.ZipFile(archive) as z:
 metadata=z.read('metadata.txt').decode()
 with z.open('results.tsv') as f:
  reader=csv.DictReader(io.TextIOWrapper(f),delimiter='\t')
  rows=[]
  for row in reader:
   if row['rsid'] in ids:
    row['prototype_lists']=';'.join(ids[row['rsid']]); rows.append(row)
fields=['rsid','variant_id','a1','a2','beta','std_error','t_stat','neg_log_p_value','sample_size','prototype_lists']
with open('outputs/webgwas_sleep_fatigue_candidate_screen.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)
print(json.dumps({'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'unique_candidates':len(ids),'matched_rows':len(rows),'matched_ids':len(set(r['rsid'] for r in rows)),'nominal_p_below_005':sum(float(r['neg_log_p_value'])> -math.log10(.05) for r in rows),'metadata':metadata},indent=2))

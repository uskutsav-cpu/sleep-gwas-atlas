import csv,hashlib,json,zipfile,io
from pathlib import Path
p=Path('work/resilient/Demographics.csv'); z=Path('work/resilient/metadata.zip')
assert hashlib.md5(p.read_bytes()).hexdigest()=='c41834a6a415afbd7422b08cb4ccfbe1'
assert hashlib.md5(z.read_bytes()).hexdigest()=='07fe55b326a6e7acc2062386d5925fc0'
with p.open(encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f))
keys=[f'phq{i}' for i in range(1,10)];complete=0;bad=0;matched=0;sum_avail=0
for row in rows:
 vals=[]
 for k in keys:
  try:v=float(row[k]); vals.append(v)
  except (ValueError,TypeError):pass
 if len(vals)==9:
  complete+=1
  if any(v not in (0,1,2,3) for v in vals):bad+=1
  try:
   total=float(row['phq_total']);sum_avail+=1;matched+=abs(sum(vals)-total)<1e-9
  except (ValueError,TypeError):pass
with zipfile.ZipFile(z) as f:
 d=list(csv.DictReader(io.StringIO(f.read('Resilient_metadata/Demographics_dictionary.csv').decode('utf-8-sig'))))
opts={r['Variable Names']:r['Categories'] for r in d if r['Variable Names'] in keys}
out={'source':'https://zenodo.org/records/16755408','data_md5':hashlib.md5(p.read_bytes()).hexdigest(),'dictionary_md5':hashlib.md5(z.read_bytes()).hexdigest(),'rows':len(rows),'nonmissing_sleep_item':sum(bool(x['phq3']) for x in rows),'nonmissing_fatigue_item':sum(bool(x['phq4']) for x in rows),'complete_9_item_rows':complete,'complete_rows_with_score_outside_0_to_3':bad,'complete_rows_with_total':sum_avail,'complete_rows_total_equal_item_sum':matched,'dictionary_phq_category_text':opts,'no_genotype_fields_in_released_dictionary':not any(any(y in ' '.join(r.values()).lower() for y in ['genotype','genomic','genetic','frailty']) for r in d),'no_participant_identifiers_in_output':True}
Path('outputs/stage6_resilient_phenotype_metadata_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='dictionary_phq_category_text'},indent=2))

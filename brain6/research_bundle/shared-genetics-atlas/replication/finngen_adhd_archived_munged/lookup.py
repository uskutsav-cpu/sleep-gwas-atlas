import csv,gzip,hashlib,json,math
from pathlib import Path
here=Path(__file__).resolve().parent
source=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged/track_b_finngen_r13_F5_ADHD.sumstats.gz')
leadpath=here.parent/'six_lead_joint/six_leads.tsv'
leads=list(csv.DictReader(leadpath.open(),delimiter='\t'));assert len(leads)==6 and len({x['variant'] for x in leads})==6
wanted={x['variant'] for x in leads}|{'rs62264764'};hits={x:[] for x in wanted};h=hashlib.sha256()
with source.open('rb') as f:
 for b in iter(lambda:f.read(1048576),b''):h.update(b)
n=0;bad=0
with gzip.open(source,'rt') as f:
 r=csv.DictReader(f,delimiter='\t');assert r.fieldnames==['SNP','A1','A2','Z','N']
 for row in r:
  n+=1
  if None in row.values() or None in row:bad+=1
  if row['SNP'] in hits:hits[row['SNP']].append(row)
results=[]
for lead in leads:
 snp=lead['variant']; rows=hits[snp];out={'variant':snp,'expected_effect':lead['effect_allele'],'other':lead['other_allele'],'matches':len(rows),'status':'MISSING' if not rows else 'DUPLICATE' if len(rows)>1 else 'ALLELES_MISMATCH'}
 if len(rows)==1:
  row=rows[0];a=row['A1'];b=row['A2'];z=float(row['Z']);p=math.erfc(abs(z)/math.sqrt(2));out.update({'file_A1':a,'file_A2':b,'file_Z':z,'N':row['N'],'p_two_sided':p})
  if {a,b}=={out['expected_effect'],out['other']} and {a,b} not in ({'A','T'},{'C','G'}):
   out['status']='MATCHED_ALLELES_UNVERIFIED_COORDINATE';out['aligned_Z']=z if a==out['expected_effect'] else -z;out['p_times_six']=min(1,6*p)
 results.append(out)
extra=hits['rs62264764']
with (here/'results.tsv').open('w',newline='') as f:
 keys=['variant','expected_effect','other','matches','status','file_A1','file_A2','file_Z','N','p_two_sided','aligned_Z','p_times_six'];w=csv.DictWriter(f,fieldnames=keys,delimiter='\t');w.writeheader();w.writerows(results)
receipt={'source':str(source),'sha256':h.hexdigest(),'rows':n,'malformed_rows':bad,'six_leads_sha256':hashlib.sha256(leadpath.read_bytes()).hexdigest(),'rs62264764_raw_matches':extra,'source_identity':'filename and historical lock only; raw source and ingest receipt not verified'}
(here/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'results':results,'receipt':receipt},indent=2))

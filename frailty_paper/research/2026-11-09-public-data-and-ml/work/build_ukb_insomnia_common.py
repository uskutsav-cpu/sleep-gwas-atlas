import csv,gzip,pathlib,hashlib,json
root=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/common_support')
dest=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity'); count=0
with gzip.open(root/'S.sumstats.gz','rt') as source,gzip.open(dest/'UKBFI.sumstats.gz','rt') as subset,gzip.open(dest/'S.sumstats.gz','wt') as out:
 sr=csv.DictReader(source,delimiter='\t');tr=csv.DictReader(subset,delimiter='\t');w=csv.DictWriter(out,fieldnames=['SNP','A1','A2','N','Z'],delimiter='\t',lineterminator='\n');w.writeheader(); target=next(tr,None)
 for row in sr:
  if target is None:break
  if row['SNP']!=target['SNP']:continue
  assert (row['A1'],row['A2'])==(target['A1'],target['A2']);w.writerow(row);count+=1;target=next(tr,None)
 assert target is None and count==1134953
print('R13 matched',count);print('S sha256',hashlib.sha256((dest/'S.sumstats.gz').read_bytes()).hexdigest())

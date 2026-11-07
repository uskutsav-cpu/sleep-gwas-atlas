import csv,gzip,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent;r=json.loads((p/'prepared_receipt.json').read_text());comp=str.maketrans('ACGT','TGCA');reference=None; checks={}
for name,digest in r['prepared_sha256'].items():
 path=p/(name+'.sumstats.gz');h=hashlib.sha256(path.read_bytes()).hexdigest();assert h==digest
 ids=hashlib.sha256();pairs=hashlib.sha256();n=0;last=''
 with gzip.open(path,'rt',newline='') as f:
  for x in csv.DictReader(f,delimiter='\t'):
   rs=x['SNP'];a=x['A1'];b=x['A2'];assert rs>last;last=rs
   assert math.isfinite(float(x['Z'])) and math.isfinite(float(x['N'])) and float(x['N'])>0
   key=min(''.join(sorted((a,b))),''.join(sorted((a.translate(comp),b.translate(comp)))))
   ids.update((rs+'\n').encode());pairs.update((rs+'\t'+key+'\n').encode());n+=1
 assert n==r['common_n'] and ids.hexdigest()==r['common_ids_sha256'];signature=pairs.hexdigest()
 if reference is None:reference=signature
 assert signature==reference
 checks[name]={'rows':n,'ordered_ids_sha256':ids.hexdigest(),'ordered_pairs_sha256':signature,'compressed_sha256':h}
(p/'input_verification.json').write_text(json.dumps({'files':checks,'result':'all files read through gzip completion with identical ordered SNP and canonical allele pair signatures'},indent=2)+'\n');print('Verified 10 complete files, identical ordered markers and alleles:',r['common_n'])

import csv,gzip,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
paths={'falls':Path('/Volumes/Brain6NHANES/neale-falls-20261028/falls.ldsc.bgz'),'questionnaire':Path('/Volumes/Brain6Work/brain6-frailty-paired/insomnia.sumstats.gz'),'registry':Path('/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity/S.sumstats.gz')}
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as x:
  for v in iter(lambda:x.read(1024*1024),b''):h.update(v)
 return h.hexdigest()
def read(f):
 d={};N=set();bad=dups=0;rows=0
 with gzip.open(f,'rt') as x:
  for row in csv.DictReader(x,delimiter='\t'):
   rows+=1;a=row['A1'];b=row['A2'];s=row['SNP'];z=float(row['Z']);n=float(row['N']);N.add(n)
   if s in d:dups+=1
   d[s]=(a,b)
   if len(a)!=1 or len(b)!=1 or a not in 'ACGT' or b not in 'ACGT' or a==b or not math.isfinite(z) or not 0<n<1e9:bad+=1
 return d,{'rows':rows,'duplicate_ids':dups,'invalid':bad,'N':sorted(N),'sha256':sha(f),'path':str(f)}
falls,first=read(paths['falls']);results={'falls':first};tr=str.maketrans('ACGT','TGCA')
for name in ('questionnaire','registry'):
 other,info=read(paths[name]);shared=falls.keys()&other.keys();same=swap=complement=pal=bad=0
 for s in shared:
  a,b=falls[s];c,d=other[s]
  if {a,b} in ({'A','T'},{'C','G'}):pal+=1
  elif (a,b)==(c,d):same+=1
  elif (a,b)==(d,c):swap+=1
  elif {a,b}=={c.translate(tr),d.translate(tr)}:complement+=1
  else:bad+=1
 info.update(shared=len(shared),direct_same=same,direct_swapped=swap,strand_compatible=complement,palindromic=pal,nonmatching=bad);results[name]=info
(p/'input_audit.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results,indent=2))

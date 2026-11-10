import csv,gzip,pathlib,json,hashlib
base=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/mak2023_same_cohort/fourway');sp=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity/S.sumstats.gz');dest=base.parent/'fiveway_R13';dest.mkdir(exist_ok=True)
sleep={};dups=set()
with gzip.open(sp,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['SNP'] in sleep:dups.add(r['SNP'])
  else:sleep[r['SNP']]=r
names=('old_FI','old_FFS','mak_FI','mak_FP');inputs=[gzip.open(base/(n+'.sumstats.gz'),'rt') for n in names];reads=[csv.DictReader(f,delimiter='\t') for f in inputs];outputs=[gzip.open(dest/(n+'.sumstats.gz'),'wt') for n in (*names,'R13')]
for f in outputs:f.write('SNP\tA1\tA2\tN\tZ\n')
comp=str.maketrans('ACGT','TGCA');n=bad=0
for rows in zip(*reads):
 first=rows[0];rs=first['SNP'];a,b=first['A1'],first['A2'];assert all(r['SNP']==rs and (r['A1'],r['A2'])==(a,b) for r in rows)
 z=sleep.get(rs)
 if z is None or rs in dups:continue
 s,t=z['A1'],z['A2']
 if (s,t)==(a,b) or (s.translate(comp),t.translate(comp))==(a,b):sgn=1
 elif (t,s)==(a,b) or (t.translate(comp),s.translate(comp))==(a,b):sgn=-1
 else:bad+=1;continue
 for f,r in zip(outputs,rows):f.write(f"{rs}\t{a}\t{b}\t{r['N']}\t{r['Z']}\n")
 outputs[-1].write(f"{rs}\t{a}\t{b}\t{z['N']}\t{float(z['Z'])*sgn:.10g}\n");n+=1
for f in inputs+outputs:f.close()
assert n>500000
meta={'n':n,'allele_mismatches_excluded':bad,'R13_duplicate_ids_excluded':len(dups),'input':'four frailty GWAS intersected with project FinnGen R13 insomnia on prior HFRS and UKB-FI restricted support, nonpalindromic, aligned alleles','sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in dest.glob('*.sumstats.gz') if not p.name.startswith('._')}}
(dest/'support.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))

import csv,gzip,hashlib,json,math
from pathlib import Path
P=Path(__file__).resolve().parent; out=Path('/Volumes/Brain6Work/brain6-frailty-paired');out.mkdir(exist_ok=True)
s=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged/insomnia.sumstats.gz')
base=Path('/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity'); files={'insomnia':s,'FI':base/'FI.sumstats.gz','FFS':base/'FFS.sumstats.gz'}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
prior=json.loads((P.parent/'receipt.json').read_text())['source_sha256'];old=json.loads((P.parents[1]/'insomnia_frailty_same_support/receipt.json').read_text())['source_sha256']
assert sha(s)==prior[str(s)]
for t in ('FI','FFS'): assert sha(files[t])==old[str(files[t])]
comp=str.maketrans('ACGT','TGCA'); valid=set('ACGT'); ambig=[{'A','T'},{'C','G'}];
def ok(r):
 try:return all(r.get(k) for k in ('SNP','A1','A2','Z','N')) and r['A1'] in valid and r['A2'] in valid and r['A1']!=r['A2'] and {r['A1'],r['A2']} not in ambig and math.isfinite(float(r['Z'])) and math.isfinite(float(r['N'])) and float(r['N'])>0
 except ValueError:return False
S={};dups=set(); stats={'insomnia_rows':0,'duplicate_insomnia_ids':0,'pairs_checked':0,'matched':0,'no_insomnia':0,'invalid_or_incompatible':0,'direct':0,'reverse':0,'strand':0,'strand_reverse':0}
with gzip.open(s,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  stats['insomnia_rows']+=1; sid=r.get('SNP');
  if sid in S or sid in dups:
   dups.add(sid);S.pop(sid,None);continue
  if ok(r):S[sid]=(r['A1'],r['A2'],r['Z'],r['N'])
stats['duplicate_insomnia_ids']=len(dups)
handles={k:gzip.open(out/(k+'.sumstats.gz'),'wt') for k in ('insomnia','FI','FFS')}
writers={k:csv.writer(f,delimiter='\t',lineterminator='\n') for k,f in handles.items()}
for w in writers.values():w.writerow(['SNP','A1','A2','Z','N'])
seen=set()
with gzip.open(files['FI'],'rt') as a,gzip.open(files['FFS'],'rt') as b:
 for fi,ffs in zip(csv.DictReader(a,delimiter='\t'),csv.DictReader(b,delimiter='\t'),strict=True):
  stats['pairs_checked']+=1; assert all(fi[k]==ffs[k] for k in ('SNP','A1','A2'))
  sid=fi['SNP']; assert sid not in seen;seen.add(sid)
  if sid not in S:stats['no_insomnia']+=1;continue
  if not ok(fi) or not ok(ffs):stats['invalid_or_incompatible']+=1;continue
  a1,a2,z,n=S[sid]; target=(fi['A1'],fi['A2']); observed=(a1,a2)
  if observed==target: sign=1; mode='direct'
  elif observed[::-1]==target:sign=-1;mode='reverse'
  elif (a1.translate(comp),a2.translate(comp))==target:sign=1;mode='strand'
  elif (a2.translate(comp),a1.translate(comp))==target:sign=-1;mode='strand_reverse'
  else:stats['invalid_or_incompatible']+=1;continue
  stats[mode]+=1;stats['matched']+=1
  writers['insomnia'].writerow([sid,*target,format(sign*float(z),'.12g'),n])
  for name,r in [('FI',fi),('FFS',ffs)]:writers[name].writerow([sid,*target,r['Z'],r['N']])
for f in handles.values():f.close()
assert stats['matched']>500000
# Re-read all streams and check ordered keys and values.
n=0
with gzip.open(out/'insomnia.sumstats.gz','rt') as a,gzip.open(out/'FI.sumstats.gz','rt') as b,gzip.open(out/'FFS.sumstats.gz','rt') as c:
 for triple in zip(*(csv.DictReader(f,delimiter='\t') for f in (a,b,c)),strict=True):
  assert all(tuple(r[k] for k in ('SNP','A1','A2'))==tuple(triple[0][k] for k in ('SNP','A1','A2')) and ok(r) for r in triple);n+=1
assert n==stats['matched']
receipt={'stats':stats,'sources_sha256':{str(v):sha(v) for v in files.values()},'output_sha256':{str(out/(k+'.sumstats.gz')):sha(out/(k+'.sumstats.gz')) for k in handles},'ordered_keys_check':True,'caveat':'FI/FFS support was previously selected using R13 and HFRS and is not the full native support.'}
(P/'prepare_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(stats,indent=2))

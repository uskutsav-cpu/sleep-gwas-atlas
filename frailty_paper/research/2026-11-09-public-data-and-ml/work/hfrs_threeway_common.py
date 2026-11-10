import gzip,csv,math,hashlib,json, pathlib
root=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation')
out=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/hfrs_fi_ffs_common_support')
source={k:root/'common_support'/f'{k}.sumstats.gz' for k in ['FI','FFS']};source['HFRS']=root/'munged/HFRS.sumstats.gz'
comp=str.maketrans('ACGT','TGCA')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def canonical(a,b,z,ta,tb):
 if any(len(x)!=1 or x not in 'ACGT' for x in [a,b,ta,tb]) or a==b or ''.join(sorted(a+b)) in ('AT','CG'):return None
 if (a,b)==(ta,tb) or (a.translate(comp),b.translate(comp))==(ta,tb):return z
 if (b,a)==(ta,tb) or (b.translate(comp),a.translate(comp))==(ta,tb):return -z
 return None
H={};duplicates=set()
with gzip.open(source['HFRS'],'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  try:z=float(r['Z']);n=float(r['N']);assert math.isfinite(z) and math.isfinite(n) and n>0
  except (ValueError,AssertionError):continue
  if r['SNP'] in H:duplicates.add(r['SNP']);continue
  H[r['SNP']]=(r['A1'],r['A2'],z,n)
print('HFRS valid',len(H),'duplicates',len(duplicates),flush=True)
files=[gzip.open(out/(k+'.sumstats.gz'),'wt') for k in ['FI','FFS','HFRS']];wr=[csv.writer(f,delimiter='\t',lineterminator='\n') for f in files]
for w in wr:w.writerow(['SNP','A1','A2','N','Z'])
count=0;skips={}
with gzip.open(source['FI'],'rt') as f1,gzip.open(source['FFS'],'rt') as f2:
 a=csv.DictReader(f1,delimiter='\t');b=csv.DictReader(f2,delimiter='\t')
 for x,y in zip(a,b):
  assert x['SNP']==y['SNP'] and (x['A1'],x['A2'])==(y['A1'],y['A2'])
  rs=x['SNP'];z=H.get(rs)
  if not z or rs in duplicates:continue
  z2=canonical(z[0],z[1],z[2],x['A1'],x['A2'])
  if z2 is None:continue
  for r in (x,y):assert math.isfinite(float(r['Z'])) and float(r['N'])>0
  wr[0].writerow([rs,x['A1'],x['A2'],x['N'],x['Z']]);wr[1].writerow([rs,y['A1'],y['A2'],y['N'],y['Z']]);wr[2].writerow([rs,x['A1'],x['A2'],z[3],z2]);count+=1
for f in files:f.close()
assert count>=200000
meta={'description':'FI and FFS originally shared insomnia-R13 support further intersected with HFRS; valid aligned nonpalindromic SNPs','n_common_snps':count,'input_sha256':{k:sha(p) for k,p in source.items()},'output_sha256':{k:sha(out/(k+'.sumstats.gz')) for k in ['FI','FFS','HFRS']},'n_hfrs_duplicate_rsIDs_excluded':len(duplicates)}
(out/'support.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta,indent=2),flush=True)

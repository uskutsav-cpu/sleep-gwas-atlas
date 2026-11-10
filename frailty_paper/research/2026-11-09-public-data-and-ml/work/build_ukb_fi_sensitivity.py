import csv,gzip,pathlib,math,hashlib,json
prev=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/hfrs_fi_ffs_common_support');out=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity')
source=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/v9-public-data/atkins-2021/Atkins_2021_FrailtyIndex_GWAS_EuropeansAged60to70years_GRCh37_UKBiobank_TwinGene_associations.txt.gz')
ref=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/eur_w_ld_chr');comp=str.maketrans('ACGT','TGCA')
target={}
with gzip.open(prev/'FI.sumstats.gz','rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):target[r['SNP']]=(r['A1'],r['A2'])
coords={}
for chr in range(1,23):
 with gzip.open(ref/f'{chr}.l2.ldscore.gz','rt') as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r['SNP'] in target:coords[r['SNP']]=(r['CHR'],r['BP'])
print('target',len(target),'coords',len(coords),flush=True)
found={};dup=set();n=0;unmatched=0;coordmiss=0
with gzip.open(source,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  n+=1;rs=r['SNP']; t=target.get(rs)
  if not t:continue
  if coords.get(rs)!=(r['CHR'],r['BP']):coordmiss+=1;continue
  try:b=float(r['BETA_UKB1']);se=float(r['SE_UKB1']);assert math.isfinite(b) and math.isfinite(se) and se>0
  except (ValueError,AssertionError):continue
  a,o=r['EA'],r['OA'];a1,a2=t
  if (a,o)==(a1,a2) or (a.translate(comp),o.translate(comp))==(a1,a2):z=b/se
  elif (o,a)==(a1,a2) or (o.translate(comp),a.translate(comp))==(a1,a2):z=-b/se
  else:unmatched+=1;continue
  if rs in found:dup.add(rs);continue
  found[rs]=z
print('read',n,'found',len(found),'duplicates',len(dup),'coordinate_mismatch',coordmiss,'allele_mismatch',unmatched,flush=True)
for rs in dup:found.pop(rs,None)
keys=['FI','FFS','HFRS'];inputs=[gzip.open(prev/(k+'.sumstats.gz'),'rt') for k in keys];outs=[gzip.open(out/(k+'.sumstats.gz'),'wt') for k in keys+['UKBFI']];readers=[csv.DictReader(f,delimiter='\t') for f in inputs];writers=[csv.writer(f,delimiter='\t',lineterminator='\n') for f in outs]
for w in writers:w.writerow(['SNP','A1','A2','N','Z'])
count=0
for rows in zip(*readers):
 x=rows[0];rs=x['SNP'];assert all(r['SNP']==rs and (r['A1'],r['A2'])==(x['A1'],x['A2']) for r in rows)
 if rs not in found:continue
 for i,r in enumerate(rows):writers[i].writerow([rs,r['A1'],r['A2'],r['N'],r['Z']])
 writers[3].writerow([rs,x['A1'],x['A2'],164610,found[rs]]);count+=1
for f in inputs+outs:f.close()
print('common',count,flush=True)
assert count>=200000
meta={'common_snps':count,'source_rows':n,'source_cohort':'Atkins UKB1, effect on reported EA, N=164610; actual per variant N unverified','matching':'rsID, GRCh37 chromosome and BP with LD reference, matching allele pair after complement/swap; duplicate source rsIDs excluded','raw_source_file':str(source),'n_coord_mismatches':coordmiss,'n_allele_mismatches':unmatched,'n_duplicates':len(dup),'output_sha256':{k:hashlib.sha256((out/(k+'.sumstats.gz')).read_bytes()).hexdigest() for k in keys+['UKBFI']}}
(out/'support.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta,indent=2),flush=True)

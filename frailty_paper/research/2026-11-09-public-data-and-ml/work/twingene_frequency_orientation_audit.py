import gzip,hashlib,json,math,collections,pathlib,statistics
src=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/v9-public-data/atkins-2021/Atkins_2021_FrailtyIndex_GWAS_EuropeansAged60to70years_GRCh37_UKBiobank_TwinGene_associations.txt.gz')
sample=lambda s:int.from_bytes(hashlib.blake2b(s.encode(),digest_size=8,person=b'FIadapter-v1').digest(),'big')%64==0
records={};dups=set();counts=collections.Counter()
with gzip.open(src,'rt') as f:
 header=next(f).strip().split('\t'); ix={name:header.index(name) for name in ('SNP','EA','OA','EAF_UKB','EAF_TG')}
 for line in f:
  counts['total_rows']+=1
  a=line.rstrip('\n').split('\t');rs=a[ix['SNP']]
  if not sample(rs):continue
  counts['sampled_rows']+=1
  if rs in records:dups.add(rs);continue
  try:
   ea,oa=a[ix['EA']].upper(),a[ix['OA']].upper();u=float(a[ix['EAF_UKB']]);t=float(a[ix['EAF_TG']]);assert ea in 'ACGT' and oa in 'ACGT' and len(ea)==1 and len(oa)==1 and ea!=oa and {ea,oa} not in ({'A','T'},{'C','G'}) and all(math.isfinite(x) and 0<=x<=1 for x in (u,t))
  except (ValueError,AssertionError,IndexError):counts['invalid_or_ambiguous']+=1;continue
  records[rs]=(u,t)
for rs in dups:records.pop(rs,None)
counts['duplicate_sampled_rsid_excluded']=len(dups);counts['unique_valid']=len(records)
def group(vals):
 d=[abs(u-t) for u,t in vals];c=[abs(u-(1-t)) for u,t in vals];n=len(vals)
 return {'n':n,'mean_abs_diff_reported':sum(d)/n,'mean_abs_diff_flipped':sum(c)/n,'median_abs_diff_reported':statistics.median(d),'median_abs_diff_flipped':statistics.median(c),'reported_closer':sum(a<b for a,b in zip(d,c)),'flipped_closer':sum(a>b for a,b in zip(d,c)),'equal':sum(a==b for a,b in zip(d,c))}
values=list(records.values());informative=[(u,t) for u,t in values if u<=.4 or u>=.6]
res={'counts':dict(counts),'all_valid':group(values),'ukb_eaf_outside_0.4_0.6':group(informative),'method':'Deterministic rsID-only hash 1/64; source cohort file non-palindromic SNPs with finite 0..1 frequencies; exclude eligible duplicate sampled rsIDs; compare |EAF_UKB-EAF_TG| to |EAF_UKB-(1-EAF_TG)|. Descriptive, no causal or probabilistic inference.','source_sha256_from_existing_archive_audit':'c7ddc6e463a212c4d10ccd938970e83c9f18ecf88b4882c56f686f272a757895','limitations':'UKB and TwinGene are different populations with differing frequencies; rare variants, sample QC and selection influence results. Frequencies may themselves be processed or inverted and do not independently prove beta sign or effect allele. LD means markers not independent.'}
pathlib.Path('outputs/twingene_frequency_orientation_audit.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))

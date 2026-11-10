import csv,gzip,hashlib,json,math,collections,pathlib
meta='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/harmonized/frailty.harmonized.tsv.gz'
cohort='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/v9-public-data/atkins-2021/Atkins_2021_FrailtyIndex_GWAS_EuropeansAged60to70years_GRCh37_UKBiobank_TwinGene_associations.txt.gz'
# sampling depends only on variant ID, not p, beta, allele or outcome
sample=lambda s: int.from_bytes(hashlib.blake2b(s.encode(),digest_size=8,person=b'FIadapter-v1').digest(),'big')%64==0
m={};dupm=set();stats=collections.Counter()
with gzip.open(meta,'rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  stats['meta_rows']+=1;rs=x['SNP']
  if not sample(rs):continue
  stats['meta_hash_sampled']+=1
  if rs in m:dupm.add(rs)
  else:m[rs]=x
for s in dupm:m.pop(s,None)
stats['meta_duplicate_sampled_rsid_excluded']=len(dupm)
records={};dupc=set()
with gzip.open(cohort,'rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  stats['cohort_rows']+=1;rs=x['SNP']
  if rs not in m:continue
  y=m[rs]
  if x['CHR']!=y['CHR'] or x['BP']!=y['BP']:stats['coordinate_mismatches']+=1;continue
  if {x['EA'],x['OA']}!={y['A1'],y['A2']} or x['EA']==x['OA']:stats['allele_mismatches']+=1;continue
  if rs in records:dupc.add(rs);continue
  try:
   b1,s1,b2,s2,bm=map(float,(x['BETA_UKB1'],x['SE_UKB1'],x['BETA_TG'],x['SE_TG'],y['BETA']))
   if not all(map(math.isfinite,(b1,s1,b2,s2,bm))) or s1<=0 or s2<=0:raise ValueError()
  except ValueError:stats['invalid_or_missing_beta_se']+=1;continue
  sign=1 if x['EA']==y['A1'] else -1; u=b1*sign;t=b2*sign;w1=1/s1**2;w2=1/s2**2
  plus=(u*w1+t*w2)/(w1+w2);minus=(u*w1-t*w2)/(w1+w2)
  records[rs]=(bm,plus,minus,abs(plus-minus),x['CHR'])
for rs in dupc:records.pop(rs,None)
stats['cohort_duplicate_eligible_rsid_excluded']=len(dupc);stats['valid_unique_matches']=len(records)
result={'method':'Deterministic BLAKE2b rsID hash modulo 64 sampled without using effect or outcome, matching original cohort file to previously harmonized FI meta file by rsID, GRCh37 chr:bp and unordered allele strings; compare inverse-variance predictions of meta beta on A1 using plus versus reversed TwinGene sign. No author confirmation.','counts':dict(stats),'meta_path':meta,'cohort_path':cohort,'groups':{}}
for threshold in (0,0.0001,0.0002,0.001):
 vals=[x for x in records.values() if x[3]>=threshold]; n=len(vals); ep=[abs(b-p) for b,p,q,d,c in vals]; em=[abs(b-q) for b,p,q,d,c in vals]
 result['groups'][str(threshold)]={'n':n,'minus_lower_error':sum(v<u for u,v in zip(ep,em)),'plus_lower_error':sum(u<v for u,v in zip(ep,em)),'ties':sum(u==v for u,v in zip(ep,em)),'mean_abs_error_plus':sum(ep)/n if n else None,'mean_abs_error_reversed':sum(em)/n if n else None,'within_0.0001_plus':sum(x<=0.0001 for x in ep),'within_0.0001_reversed':sum(x<=0.0001 for x in em)}
print(json.dumps(result,indent=2));pathlib.Path('outputs/twingene_direction_genomewide_audit.json').write_text(json.dumps(result,indent=2)+'\n')

import csv,gzip,math,pathlib,collections,json,hashlib,itertools
root=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/mak2023_same_cohort'); ref=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/hm3_grch37_variant_map.tsv.gz'); stats=collections.Counter(); mapping={};duplicates=set()
with gzip.open(ref,'rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  if x['SNP'] in mapping:duplicates.add(x['SNP'])
  else:mapping[x['SNP']]=(x['CHR'],x['BP'],x['A1'],x['A2'])
for k in duplicates:mapping.pop(k,None)
results={};dups=set()
with (root/'fi.fastGWA').open() as f1,(root/'fp.fastGWA').open() as f2:
 for a,b in itertools.zip_longest(csv.DictReader(f1,delimiter='\t'),csv.DictReader(f2,delimiter='\t')):
  stats['source_rows']+=1
  if a is None or b is None:stats['length_mismatch']+=1;continue
  if (a['SNP'],a['CHR'],a['POS'])!=(b['SNP'],b['CHR'],b['POS']):stats['row_order_mismatch']+=1;continue
  rs=a['SNP'];m=mapping.get(rs)
  if m is None:continue
  if rs in results:dups.add(rs);continue
  chrom,pos,ra,rb=m
  if a['CHR']!=chrom or a['POS']!=pos:stats['reference_position_mismatch']+=1;continue
  A,B=a['A1'],a['A2'];AA,BB=b['A1'],b['A2']
  if len(A)!=1 or len(B)!=1 or A not in 'ACGT' or B not in 'ACGT' or A==B or {A,B} in ({'A','T'},{'C','G'}):stats['ambiguous_or_invalid']+=1;continue
  if {A,B}!={ra,rb} or {A,B}!={AA,BB}:stats['allele_mismatch']+=1;continue
  try:
   n1,n2=float(a['N']),float(b['N']); beta1,se1,beta2,se2=map(float,(a['BETA'],a['SE'],b['BETA'],b['SE']));af1,af2=float(a['AF1']),float(b['AF1'])
   if not all(math.isfinite(v) for v in (n1,n2,beta1,se1,beta2,se2,af1,af2)) or min(se1,se2)<=0 or min(n1,n2)<403041*2/3 or not (0<af1<1 and 0<af2<1):raise ValueError()
  except (ValueError,TypeError):stats['numeric_qc_failed']+=1;continue
  sign=1 if (AA,BB)==(A,B) else -1
  results[rs]=(A,B,n1,beta1/se1,n2,sign*beta2/se2)
for rs in dups:results.pop(rs,None)
stats['duplicate_eligible_ids_excluded']=len(dups);stats['retained']=len(results)
assert stats['length_mismatch']==0 and stats['row_order_mismatch']==0
assert len(results)>300000
with gzip.open(root/'mak_fi_common.sumstats.gz','wt') as fi,gzip.open(root/'mak_fp_common.sumstats.gz','wt') as fp:
 fi.write('SNP\tA1\tA2\tN\tZ\n');fp.write('SNP\tA1\tA2\tN\tZ\n')
 for rs in sorted(results):
  a,b,n1,z1,n2,z2=results[rs];fi.write(f'{rs}\t{a}\t{b}\t{n1:.0f}\t{z1:.10g}\n');fp.write(f'{rs}\t{a}\t{b}\t{n2:.0f}\t{z2:.10g}\n')
summary={'counts':dict(stats),'notes':'Exploratory common support; exact nonpalindromic GRCh37 position+allele match, unique IDs, per-variant N >= 2/3 of 403041, Z=beta/se, no phenotype scores; reference map construction and source build still under review.','source_md5_verified':{'fi.fastGWA':'5cd10928bde27e48b82d38f097696473','fp.fastGWA':'3ec3de20038f3022eaa9163347b0e194'},'allele_frequency_check':'0<AF1<1; no AF difference filter between sources; same input alleles verified; reference ancestry and AF mismatch not separately assessed.'}
pathlib.Path('outputs/mak2023_same_cohort_preparation.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

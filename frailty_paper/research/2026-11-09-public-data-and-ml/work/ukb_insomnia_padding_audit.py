import csv,gzip,collections,json,pathlib,math
src=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/munged/insomnia.sumstats.gz');common=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/mak2023_same_cohort/fourway/mak_FI.sumstats.gz');stats=collections.Counter();blank=set();observed=set();examples=[]
with gzip.open(src,'rt') as f:
 for row in csv.DictReader(f,delimiter='\t'):
  stats['rows']+=1;is_blank=all(not row[x] for x in ('A1','A2','Z','N'))
  if is_blank:
   stats['all_four_fields_blank']+=1;blank.add(row['SNP'])
   if len(examples)<5:examples.append(row['SNP'])
  elif all(row[x] for x in ('A1','A2','Z','N')):stats['all_four_fields_present']+=1;observed.add(row['SNP'])
  else:stats['other_missingness']+=1
seen=collections.Counter()
with gzip.open(common,'rt') as f:
 for row in csv.DictReader(f,delimiter='\t'):
  rs=row['SNP'];seen['common']+=1
  seen['all_four_fields_blank' if rs in blank else 'all_four_fields_present' if rs in observed else 'absent_or_partial']+=1
res={'all_munged':dict(stats),'fourway_overlap':dict(seen),'example_blank_ids':examples,'origin':'Existing munge log reports writing the full allele-merge list, including rows without beta. Fully blank statistic rows are consistent with such padding, not demonstrated raw-data allele defects.','note':'Do not fill rows with inferred or reference alleles: the beta, Z and N are also missing. This check does not independently reconstruct the harmonized data or count reasons for absence.'};pathlib.Path('outputs/ukb_insomnia_padding_audit.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))

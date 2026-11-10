import csv,gzip,collections,json,pathlib
p='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/munged/insomnia.sumstats.gz';q='/Volumes/Extreme SSD/Codex Archive/mak2023_same_cohort/fourway/mak_FI.sumstats.gz'
with gzip.open(p,'rt') as f:s={x['SNP']:(x['A1'],x['A2']) for x in csv.DictReader(f,delimiter='\t')}
c=str.maketrans('ACGT','TGCA');counts=collections.Counter()
with gzip.open(q,'rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  if x['SNP'] not in s:counts['absent']+=1;continue
  a,b=x['A1'],x['A2'];t,u=s[x['SNP']]
  if not t or not u:counts['blank_allele_fields']+=1
  elif set((a,b)) in (set((t,u)),set((t.translate(c),u.translate(c)))):counts['matching_alleles']+=1
  else:counts['other_incompatible_allele_fields']+=1
out=dict(counts);pathlib.Path('outputs/mak_ukb_insomnia_allele_qc.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

import csv,gzip,pathlib,hashlib,json,collections
pgs=pathlib.Path('work/PGS005229_complete.txt.gz');src=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv')
target={};othermissing=0
with gzip.open(pgs,'rt',newline='') as f:
 for x in csv.DictReader((v for v in f if not v.startswith('#')),delimiter='\t'):
  if not x['chr_name'] or not x['chr_position']:
   if x['rsID'] in target:raise RuntimeError('duplicate score rsID')
   target[x['rsID']]=x['effect_allele'];othermissing+=not bool(x['other_allele'])
assert len(target)==388036
found=collections.defaultdict(list);n=0
with src.open(newline='') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  n+=1
  if x['SNP'] in target:found[x['SNP']].append((x['CHROM'],x['POS'],x['A1'],x['A2']))
c=collections.Counter()
for id,allele in target.items():
 x=found.get(id,[])
 if not x:c['rsid_absent_in_FFS_source']+=1
 elif len(x)>1:c['rsid_multiple_rows_unresolved']+=1
 elif allele==x[0][2]:c['rsid_unique_effect_matches_source_A1_no_score_coord']+=1
 elif allele==x[0][3]:c['rsid_unique_effect_matches_source_A2_no_score_coord']+=1
 else:c['rsid_unique_effect_not_in_source_pair']+=1
assert sum(c.values())==len(target)
result={'score_id':'PGS005229','score_rows_missing_coordinate':len(target),'of_these_missing_other_allele':othermissing,'FFS_source_rows':n,'rsid_only_counts':dict(c),'limitations':'Score coordinates absent; rsID-only matching is not a verified positional or biallelic match; no missing coordinates or alleles imputed; no score or prediction evaluated.'}
pathlib.Path('outputs/ffs_pgs_missing_coordinate_rsid_audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

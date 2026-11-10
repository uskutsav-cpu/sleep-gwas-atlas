import csv,gzip,math,pathlib,json,statistics
src='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/v9-public-data/atkins-2021/Atkins_2021_FrailtyIndex_GWAS_EuropeansAged60to70years_GRCh37_UKBiobank_TwinGene_associations.txt.gz'
rows=list(csv.DictReader(open('outputs/expanded_22_candidate_risk_gene_table.csv')));wanted={x['variant']:x for x in rows};refs={x['variant']:x for x in csv.DictReader(open('outputs/reference_22_score_variance_decomposition.csv'))};hits={}
with gzip.open(src,'rt') as f:
 header=next(f).rstrip('\n').split('\t');ix={n:header.index(n) for n in ('SNP','CHR','BP','EA','OA','EAF_UKB','EAF_TG')}
 for line in f:
  a=line.rstrip('\n').split('\t');rs=a[ix['SNP']]
  if rs not in wanted:continue
  x=wanted[rs]
  if a[ix['CHR']]!=x['GRCh37_chr_assumed'] or a[ix['BP']]!=x['GRCh37_pos_assumed'] or {a[ix['EA']],a[ix['OA']]}!={x['counted_allele'],x['other_allele']}:continue
  assert rs not in hits
  u=float(a[ix['EAF_UKB']]);t=float(a[ix['EAF_TG']]);ref=float(refs[rs]['reference_counted_allele_frequency']);assert 0<=u<=1 and 0<=t<=1
  if a[ix['EA']]!=x['counted_allele']:u=1-u;t=1-t
  hits[rs]={'variant':rs,'counted_allele':x['counted_allele'],'reference_counted_allele_frequency':ref,'UKB_frequency_assuming_source_EAF_UKB_is_source_EA':u,'TwinGene_frequency_assuming_source_EAF_TG_is_source_EA':t,'TwinGene_frequency_if_source_EAF_TG_is_source_OA':1-t,'UKB_vs_reference_absdiff':abs(u-ref),'TG_reported_vs_reference_absdiff':abs(t-ref),'TG_inverted_vs_reference_absdiff':abs(1-t-ref)}
assert len(hits)==22
out=pathlib.Path('outputs/22_candidate_cohort_allele_frequency_check.csv');dat=[hits[x['variant']] for x in rows]
with out.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=dat[0]);w.writeheader();w.writerows(dat)
report={'n':len(dat),'mean_abs_difference_UKB_to_1000G':statistics.mean(r['UKB_vs_reference_absdiff'] for r in dat),'mean_abs_difference_TG_reported_to_1000G':statistics.mean(r['TG_reported_vs_reference_absdiff'] for r in dat),'mean_abs_difference_TG_inverted_to_1000G':statistics.mean(r['TG_inverted_vs_reference_absdiff'] for r in dat),'inverted_closer_count':sum(r['TG_inverted_vs_reference_absdiff']<r['TG_reported_vs_reference_absdiff'] for r in dat),'reported_closer_count':sum(r['TG_inverted_vs_reference_absdiff']>r['TG_reported_vs_reference_absdiff'] for r in dat),'equal_count':sum(r['TG_inverted_vs_reference_absdiff']==r['TG_reported_vs_reference_absdiff'] for r in dat),'note':'Descriptive; 1000G EUR, UKB, and TG differ in ancestry, sampling and imputation; frequency not proof of beta orientation.'}
pathlib.Path('outputs/22_candidate_cohort_allele_frequency_check.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

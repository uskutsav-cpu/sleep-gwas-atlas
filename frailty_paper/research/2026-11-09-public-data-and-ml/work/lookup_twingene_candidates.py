import csv,gzip,math,json
from collections import defaultdict
src='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/v9-public-data/atkins-2021/Atkins_2021_FrailtyIndex_GWAS_EuropeansAged60to70years_GRCh37_UKBiobank_TwinGene_associations.txt.gz'
with open('outputs/unvalidated_15_variant_score_prototype.csv',newline='') as f: candidates=list(csv.DictReader(f))
ids={r['variant'] for r in candidates}; found=defaultdict(list)
with gzip.open(src,'rt') as f:
 for row in csv.DictReader(f,delimiter='\t'):
  if row['SNP'] in ids: found[row['SNP']].append(row)
fields=['variant','chr_GRCh37','pos_GRCh37','counted_allele','source_rows_for_rsid','matching_rows','source_effect_allele','source_other_allele','UKB_beta_aligned','UKB_se','UKB_p','TG_beta_raw','TG_beta_aligned_adapter','TG_se','TG_ci95_low','TG_ci95_high','TG_p','identity_and_interpretation']
with open('outputs/candidate_twingene_fi_cohort_lookup.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
 for c in candidates:
  rs=c['variant']; hits=found[rs]; same=[r for r in hits if r['CHR']==c['chr'] and r['BP']==c['position_GRCh37_inferred'] and {r['EA'],r['OA']}=={c['counted_allele'],c['other_allele']}]
  out=dict(variant=rs,chr_GRCh37=c['chr'],pos_GRCh37=c['position_GRCh37_inferred'],counted_allele=c['counted_allele'],source_rows_for_rsid=len(hits),matching_rows=len(same),identity_and_interpretation='Exploratory; source cohort not independent of meta-GWAS; alleles/representations require verification; no clinical inference')
  if len(same)==1:
   r=same[0];sign=1 if r['EA']==c['counted_allele'] else -1
   out.update(source_effect_allele=r['EA'],source_other_allele=r['OA'],UKB_p=r['P_UKB1'],TG_beta_raw=r['BETA_TG'],TG_se=r['SE_TG'],TG_p=r['P_TG'])
   try: out['UKB_beta_aligned']=sign*float(r['BETA_UKB1']);out['UKB_se']=float(r['SE_UKB1'])
   except ValueError:pass
   try:
    b=-sign*float(r['BETA_TG']);se=float(r['SE_TG']);out.update(TG_beta_aligned_adapter=b,TG_ci95_low=b-1.95996398454*se,TG_ci95_high=b+1.95996398454*se)
   except ValueError:pass
  w.writerow(out)
print(json.dumps({'candidates':len(candidates),'source_rsid_hit_counts':{k:len(v) for k,v in found.items()}}))

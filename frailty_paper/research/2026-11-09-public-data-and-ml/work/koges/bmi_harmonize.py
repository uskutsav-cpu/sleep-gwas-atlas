import csv
r={x['variant']:x for x in csv.DictReader(open('outputs/stage6_koges_bmi_candidate_followup.csv'))}
rows=[]
for fn,kind in [('outputs/unvalidated_15_variant_score_prototype.csv','15'),('outputs/unvalidated_22_shared_direction_score.csv','22-count')]:
 for x in csv.DictReader(open(fn)):
  b=r[x['variant']]
  if b['status']!='exact pair found' or float(b['BH_q_33'])>=.05:continue
  allele=x['counted_allele'];other=x['other_allele']; assert {allele,other}=={b['ref'],b['alt']}
  beta=float(b['beta_alt'])*(1 if allele==b['alt'] else -1)
  row={'score':kind,'variant':x['variant'],'counted_allele':allele,'BMI_beta_for_counted_allele_rounded':beta,'BMI_p_rounded':b['p'],'BMI_BH_q_33_from_rounded_p':b['BH_q_33']}
  if kind=='15':
   for key in ('FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight'):row[key]=x[key]
  rows.append(row)
assert len(rows)==4
fields=['score','variant','counted_allele','BMI_beta_for_counted_allele_rounded','BMI_p_rounded','BMI_BH_q_33_from_rounded_p','FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight']
with open('outputs/stage6_koges_bmi_score_allele_harmonization.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

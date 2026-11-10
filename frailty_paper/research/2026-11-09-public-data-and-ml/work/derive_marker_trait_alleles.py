import csv
src=list(csv.DictReader(open('outputs/candidate_risk_gene_review_15.csv')))
rows=[]
for r in src:
 a,b=r['counted_allele'],r['other_allele']
 bs=[float(r[k]) for k in ('FI_meta_beta','FFS_beta_inferred_A1','insomnia_Jansen_beta')]
 ps=[float(r[k]) for k in ('FI_meta_p','FFS_p','insomnia_Jansen_p')]
 alleles=[a if x>0 else b if x<0 else 'undetermined' for x in bs]
 row={'variant':r['variant'],'GRCh37_chr':r['chr_GRCh37'],'GRCh37_pos':r['position_GRCh37'],'reference_multiallelic':r['multiallelic_reference'],'nearest_gene_body':r['nearest_gene_body'],'counted_allele_in_source':a,'other_allele_in_source':b,'FI_higher_direction_allele':alleles[0],'FFS_higher_direction_allele':alleles[1],'insomnia_higher_direction_allele':alleles[2],'three_direction_alleles_same':len(set(alleles))==1,'FI_nominal_p':ps[0],'FFS_nominal_p':ps[1],'insomnia_nominal_p':ps[2],'all_three_nominal_p_below_0_05':all(p<.05 for p in ps),'FI_ukb_vs_twingene_direction':r['FI_cohort_direction_check'],'interpretation':'Estimated directional allele in reported two-allele contrast only; post-selected, unvalidated; no absolute risk; opposite-allele flipping requires biallelic genotype representation'}
 rows.append(row)
with open('outputs/candidate_trait_direction_alleles.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print('rows',len(rows),'same directions',sum(x['three_direction_alleles_same'] for x in rows),'all three nominal',sum(x['all_three_nominal_p_below_0_05'] for x in rows),'both',sum(x['all_three_nominal_p_below_0_05'] and x['three_direction_alleles_same'] for x in rows))

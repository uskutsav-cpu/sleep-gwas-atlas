import csv, math
files={'base':'outputs/candidate_risk_gene_review_15.csv','hfrs':'outputs/candidate_finngen_hfrs_full_gwas_lookup.csv','direction':'outputs/candidate_trait_direction_alleles.csv','vep':'outputs/candidate_vep_allele_consequences.csv','ld':'outputs/candidate_nearest_ffs_lead_ld.csv','five':'outputs/prototype_five_source_descriptive_triage.csv'}
load=lambda p:{r['variant']:r for r in csv.DictReader(open(p))}
base=load(files['base']);h=load(files['hfrs']);d=load(files['direction']);vep=load(files['vep']);five=load(files['five'])
ld={r['candidate']:r for r in csv.DictReader(open(files['ld']))}
assert all(set(y)==set(base) for y in [h,d,vep,five,ld]);rows=[]
for rs,b in base.items():
 x=h[rs];z=d[rs];v=vep[rs];l=ld[rs];f=five[rs]
 flags=[]
 if b['multiallelic_reference']=='yes':flags.append('reference_multiallelic')
 if b['FI_cohort_direction_check']=='opposite':flags.append('FI_component_direction_discordance')
 if z['three_direction_alleles_same']=='False':flags.append('FI_FFS_UKB_insomnia_directions_differ')
 if l['reference_LD_status']=='returned_rsID_based':flags.append('near_FFS_lead_has_reference_LD')
 if rs=='rs8067974':flags.append('published_brain_HEXIM1_QTL_lead')
 if rs=='rs12536335':flags.append('published_HFRS_without_dementia_hit')
 if rs=='rs7734045':flags.append('near_TRPC7_insomnia_gene_context')
 row={'variant':rs,'GRCh37_chr':b['chr_GRCh37'],'GRCh37_pos':b['position_GRCh37'],'GRCh38_pos':x['GRCh38_pos'],'counted_allele':b['counted_allele'],'other_allele':b['other_allele'],'nearest_gene_body_positional':b['nearest_gene_body'],'distance_bp':b['distance_to_gene_body_bp'],'FI_beta_on_counted_allele':b['FI_meta_beta'],'FI_p':b['FI_meta_p'],'FFS_beta_inferred_on_counted_allele':b['FFS_beta_inferred_A1'],'FFS_p':b['FFS_p'],'UKB_insomnia_beta_on_counted_allele':b['insomnia_Jansen_beta'],'UKB_insomnia_p':b['insomnia_Jansen_p'],'FinnGen_insomnia_beta_on_counted_allele':b['insomnia_FinnGen_beta'],'FinnGen_insomnia_p':b['insomnia_FinnGen_p'],'MVP_insomnia_beta_on_counted_allele':b['insomnia_MVP_beta'],'MVP_insomnia_p':b['insomnia_MVP_p'],'FinnGen_HFRS_beta_on_counted_allele':x['hfrs_beta_on_counted_allele'],'FinnGen_HFRS_se':x['hfrs_se'],'FinnGen_HFRS_p':x['hfrs_p'],'five_input_same_direction_nominal_p_below_0_05':str(f['all_five_nominal_p_lt_0_05']=='True' and f['all_five_signs_equal']=='True'),'HFRS_nominal_p_below_0_05':str(float(x['hfrs_p'])<.05),'VEP_transcript_terms_for_contrast_alternate':v['transcript_terms'],'closest_FFS_lead_if_any':l['nearest_ffs_lead'],'FFS_lead_reference_r2_if_returned':l['reference_r2'],'flags':';'.join(flags),'interpretation':'postselected descriptive research candidates; nonindependent inputs; no confirmed shared causal variant, score validation or individual risk'}
 rows.append(row)
with open('outputs/integrated_15_candidate_review.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
print('rows',len(rows),'five input same nominal',sum(r['five_input_same_direction_nominal_p_below_0_05']=='True' for r in rows),'HFRS nominal',sum(r['HFRS_nominal_p_below_0_05']=='True' for r in rows))

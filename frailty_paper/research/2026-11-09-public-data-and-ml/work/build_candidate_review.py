import csv
from pathlib import Path

def read(p):
 with open(p,newline='') as f:return list(csv.DictReader(f))
a=read('outputs/prototype_variant_five_trait_evidence.csv');b=read('outputs/candidate_nearest_gene_annotations.csv');c=read('outputs/candidate_twingene_fi_cohort_lookup.csv');d=read('outputs/unvalidated_15_variant_score_prototype.csv')
def index(xs):
 result={r['variant']:r for r in xs};assert len(result)==len(xs);return result
b,c,d=map(index,(b,c,d));assert len(a)==15 and set(b)==set(c)==set(d)=={x['variant'] for x in a}
fields=['variant','chr_GRCh37','position_GRCh37','counted_allele','other_allele','multiallelic_reference','cross_trait_source_rank','nearest_gene_body','nearest_gene_ensembl_ids','distance_to_gene_body_bp','FI_meta_beta','FI_meta_se','FI_meta_p','UKB_FI_beta','UKB_FI_se','UKB_FI_p','TwinGene_FI_beta_adapter','TwinGene_FI_se','TwinGene_FI_p','FI_cohort_direction_check','FFS_beta_inferred_A1','FFS_se','FFS_p','insomnia_Jansen_beta','insomnia_Jansen_se','insomnia_Jansen_p','insomnia_FinnGen_beta','insomnia_FinnGen_se','insomnia_FinnGen_p','insomnia_MVP_beta','insomnia_MVP_p','gene_annotation_source','caution']
out=[]
for q in a:
 v=q['variant'];g=b[v];t=c[v];p=d[v]; assert q['counted_allele']==g['counted_allele']==t['counted_allele']==p['counted_allele']; assert q['position_GRCh37']==g['position_GRCh37']==t['pos_GRCh37']==p['position_GRCh37_inferred'];assert t['matching_rows']=='1'
 try:
  u=float(t['UKB_beta_aligned']);w=float(t['TG_beta_aligned_adapter']);direction='same' if u*w>0 else 'opposite' if u*w<0 else 'zero'
 except ValueError:direction='unavailable'
 row=dict(variant=v,chr_GRCh37=q['chromosome_GRCh37'],position_GRCh37=q['position_GRCh37'],counted_allele=q['counted_allele'],other_allele=q['other_allele'],multiallelic_reference=q['reference_multiallelic'],cross_trait_source_rank=q['published_cross_trait_rank'],nearest_gene_body=g['nearest_gene_names'],nearest_gene_ensembl_ids=g['nearest_gene_ids'],distance_to_gene_body_bp=g['distance_to_gene_body_bp'],FI_meta_beta=q['FI_beta'],FI_meta_se=q['FI_se'],FI_meta_p=q['FI_p'],UKB_FI_beta=t['UKB_beta_aligned'],UKB_FI_se=t['UKB_se'],UKB_FI_p=t['UKB_p'],TwinGene_FI_beta_adapter=t['TG_beta_aligned_adapter'],TwinGene_FI_se=t['TG_se'],TwinGene_FI_p=t['TG_p'],FI_cohort_direction_check=direction,FFS_beta_inferred_A1=q['FFS_beta'],FFS_se=q['FFS_se'],FFS_p=q['FFS_p'],insomnia_Jansen_beta=q['Jansen_insomnia_beta'],insomnia_Jansen_se=q['Jansen_insomnia_se'],insomnia_Jansen_p=q['Jansen_insomnia_p'],insomnia_FinnGen_beta=q['R13_insomnia_beta'],insomnia_FinnGen_se=q['R13_insomnia_se'],insomnia_FinnGen_p=q['R13_insomnia_p'],insomnia_MVP_beta=q['MVP_insomnia_beta'],insomnia_MVP_p=q['MVP_insomnia_p'],gene_annotation_source=g['source_url'],caution='Exploratory post-selection; p unadjusted; cohorts overlap or overlap unknown; TG adapter inferred from audit; no proven causal gene, joint risk score, replication or clinical use')
 out.append(row)
with open('outputs/candidate_risk_gene_review_15.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
print('validated:',len(out),'rows; identical IDs, positions and counted alleles across four inputs')

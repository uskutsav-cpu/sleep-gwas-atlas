import csv,pathlib,collections,hashlib,json
p=pathlib.Path(__file__).resolve().parent;root=p.parent
read=lambda p,sep='\t':list(csv.DictReader(p.open(),delimiter=sep))
candidates=[r for r in read(root/'evidence/BRAIN6_FINAL_CANDIDATES.tsv') if r['pair_id']=='insomnia__adhd'];selected=read(root/'variant_context/selected_variant_gene_context.tsv');gt=read(root/'variant_context/gtex_v8/significant_associations.tsv');gmap=collections.defaultdict(list)
for r in gt:gmap[r['snp']].append(r)
out=[]
for c in candidates:
 x=[v for v in selected if v['region']==c['candidate_locus_id']];assert x
 choice={}
 for k in ('max_abs_adult_dnase','max_avi','max_insomnia_pip','max_adhd_pip'):
  hit=[v for v in x if k in v['selection_reasons'].split(';')];assert len(hit)==1;choice[k]=hit[0]
 brain=sorted({a['snp']+'|'+a['gene_symbol']+'|'+a['tissue'] for v in x for a in gmap[v['snp']] if a['type']=='Eqtl' and a['tissue'].startswith('Brain_')})
 y={'region_grch37':c['geographic_region_grch37'],'candidate_locus_id':c['candidate_locus_id'],'placo_lead':c['lead_variants'],'placo_lead_p':c['placo_lead_p'],'selected_union_coverage':sum(1 for a in read(root.parent/'brain6-alphagenome-pilot/results/jbrowse_expanded/expanded_variant_evidence.tsv') if a['candidate_locus_id']==c['candidate_locus_id']),'ABF_coloc_PP4_p12_1e-5':c['trait_coloc_pp_h4'],'ABF_coloc_PP4_p12_1e-6':c['trait_coloc_low_prior_pp_h4'],'lava_family_status':c['canonical_lava_family'],'alternative_local_status':c['alternative_local_status'],'qtl_coloc_status':c['eqtl_coloc_status'],'two_trait_independent_replication':c['independent_two_trait_locus_replication'],'frozen_evidence_class':c['final_evidence_class'],'gtex_v8_brain_eqtl_for_20_selected':';'.join(brain) if brain else 'NO_RETURNED_SIGNIFICANT_PAIR_IN_SELECTED_SUBSET'}
 for label,v in choice.items():
  for src in ['snp','grch38_ref_alt','insomnia_pip','adhd_pip','insomnia_credible_set_member','adhd_credible_set_member','adult_dnase_log2ratio','avi_phred','overlapping_genes']:
   y[label+'_'+src]=v[src]
 out.append(y)
assert len(out)==6
with (p/'brain6_six_region_atlas.tsv').open('w',newline='') as f:w=csv.DictWriter(f,out[0],delimiter='\t');w.writeheader();w.writerows(out)
song=read(root/'evidence/song124_mvp_hfrs_joint_screen.csv',',');genes={x['variant']:x for x in read(root/'evidence/song124_six_hit_positional_genes.csv',',')};marg={x['candidate']:x for x in read(p/'published_candidate_marginal_lookup.csv',',')};f=[]
for s in song:
 if s['eligible_for_joint_screen']!='True' or float(s['joint_bonferroni_124'])>=.05:continue
 g=genes[s['variant']];q=marg[s['variant']];f.append({'variant':s['variant'],'song_chr_assumed_grch37':s['song_chr'],'song_pos_assumed_grch37':s['song_position_reported'],'effect_allele_for_comparison':s['song_allele_assumed_38'],'other_allele':s['song_other_assumed_38'],'nearest_gene_positional_only':g['nearest_names'],'nearest_gene_id':g['nearest_ids'],'gene_distance_bp':g['distance_bp'],'reference_alleles':g['reference_alleles'],'MVP_insomnia_p':s['MVP_p'],'FinnGen_HFRS_p':s['HFRS_p'],'max_p':s['joint_p_max'],'max_p_times_124':s['joint_bonferroni_124'],'same_direction_descriptive':s['same_effect_direction'],'FFS_marginal_p':q['FFS_P'],'FFS_allele_convention':'inferred','independence':'UNKNOWN','selection':'POST_SELECTED_PUBLISHED_VARIANTS','status':'EXPLORATORY_NOT_PROVEN_SHARED_CAUSAL_VARIANT'})
assert len(f)==6
with (p/'frailty_six_candidate_atlas.tsv').open('w',newline='') as o:w=csv.DictWriter(o,f[0],delimiter='\t');w.writeheader();w.writerows(f)
print('wrote',len(out),'Brain6 regions and',len(f),'frailty candidates')

import csv,pathlib
m=list(csv.DictReader(open('outputs/shared22_eight_source_evidence_matrix.csv')));q=list(csv.DictReader(open('outputs/stage3_gtex_v10_candidate_qtls.csv')));ex=next(x for x in csv.DictReader(open('outputs/integrated_15_candidate_review.csv')) if x['variant']=='rs323509')
rows=[]
def edge(a,b,etype,source,trait,tissue,p,details):
 rows.append({'source_node':a,'target_node':b,'edge_type':etype,'evidence_source':source,'phenotype_or_molecular_trait':trait,'tissue':tissue,'p_value':p,'details':details})
for x in m:
 for key,label in [('atkins_FI','Atkins FI'),('ye_FFS','Ye FFS'),('mak_FI','Mak FI'),('mak_FP','Mak modified physical frailty'),('hfrs','FinnGen HFRS'),('ukb_insomnia','UKB insomnia'),('mvp_insomnia','MVP insomnia'),('finngen_insomnia','FinnGen insomnia')]:
  edge('variant:'+x['variant'],'phenotype:'+label,'GWAS_marginal_association',label,label,'',x[key+'_p'],'beta '+x[key+'_beta_aligned']+' counted allele '+x['score_counted_allele']+'/'+x['other_allele']+'; post-selected; not causal')
for source,trait,p,beta in [('Atkins FI','FI','FI_p','FI_beta_on_counted_allele'),('Ye FFS','FFS','FFS_p','FFS_beta_inferred_on_counted_allele'),('UKB insomnia','UKB insomnia','UKB_insomnia_p','UKB_insomnia_beta_on_counted_allele'),('MVP insomnia','MVP insomnia','MVP_insomnia_p','MVP_insomnia_beta_on_counted_allele'),('FinnGen insomnia','FinnGen insomnia','FinnGen_insomnia_p','FinnGen_insomnia_beta_on_counted_allele'),('FinnGen HFRS','HFRS','FinnGen_HFRS_p','FinnGen_HFRS_beta_on_counted_allele')]:
 edge('variant:rs323509','phenotype:'+source,'GWAS_marginal_association',source,trait,'',ex[p],'beta '+ex[beta]+' counted allele '+ex['counted_allele']+'/'+ex['other_allele']+'; source coding uncertain; post-selected; not causal')
for x in q:
 gene=(x['gencode_id'] or x['gene_symbol']).split('.')[0]; assert gene
 edge('variant:'+x['candidate'],'gene:'+gene,'bulk_tissue_'+x['qtl_type'],'GTEx portal gtex_v10',x['qtl_type'],x['tissue'],x['p_value'],'GTEx ID '+x['gtex_variant_id']+'; gene symbol '+x['gene_symbol']+'; molecular event '+x['molecular_phenotype']+'; exact-marker association not colocalization; effect direction not aligned')
for path,trait in [('outputs/stage3_smim32_finngen_coloc_abf_sensitivity.csv','FinnGen insomnia'),('outputs/stage3_smim32_mak_coloc_abf_sensitivity.csv','Mak FI'),('outputs/stage3_smim32_mak_coloc_abf_sensitivity.csv','Mak modified physical frailty')]:
 data=list(csv.DictReader(open(path))); selected=[x for x in data if float(x['p12'])==1e-5 and (('trait' not in x and float(x['qtl_prior_effect_sd'])==.15) or ('trait' in x and x['trait']==('Mak_FI' if trait=='Mak FI' else 'Mak_FP') and float(x['gwas_prior_scale_factor'])==1))];assert len(selected)==1
 x=selected[0];edge('locus:chr5:135702653-136702653_GRCh38','gene:ENSG00000271824','exploratory_region_colocalization','eQTL Catalogue r7 muscle and '+trait,trait,'skeletal muscle','','H1='+x['H1']+'; H3='+x['H3']+'; H4='+x['H4']+'; single signal model, prior-sensitive, incomplete credible set coverage; not interpretable as a definitive shared or distinct signal')
for x in csv.DictReader(open('outputs/stage3_smad5_splice_coloc_abf_sensitivity.csv')):
 if float(x['p12'])==1e-5 and float(x['qtl_prior_effect_sd'])==.15:
  edge('locus:chr5:135702653-136702653_GRCh38','gene:ENSG00000113658','exploratory_region_splice_colocalization','eQTL Catalogue r7 QTD000285 and '+x['trait'],x['trait'],'skeletal muscle','','event 5:136132962:136153592:clu_31306_+; H1='+x['H1']+'; H3='+x['H3']+'; H4='+x['H4']+'; multiple molecular signals and incomplete credible set coverage; single signal H3/H4 not definitive')
one=pathlib.Path('outputs/stage4_onek1k_candidate_credible_set_matches.csv')
if one.exists():
 seen=set()
 for x in csv.DictReader(one.open()):
  k=(x['candidate'],x['qtl_variant_id'],x['dataset'],x['molecular_trait_id'],x['credible_set_id'])
  if k in seen:continue
  seen.add(k)
  edge('variant:'+x['candidate'],'gene:'+x['gene_id'],'cell_type_eQTL_credible_set','eQTL Catalogue r7 '+x['dataset'],x['molecular_trait_id'],x['cell_population'],x['p_value'],'variant '+x['qtl_variant_id']+'; PIP '+x['pip']+'; set '+x['credible_set_id']+'; cell QTL only; no cell-specific GWAS colocalization or enrichment')
tp=pathlib.Path('outputs/stage3_published_fi_insomnia_twas_intersection.csv')
if tp.exists():
 for x in csv.DictReader(tp.open()):
  for title,key in [('FI','fi'),('insomnia','insomnia')]:
   edge('gene_symbol:'+x['gene'],'phenotype:published_'+title,'published_whole_blood_FUSION_TWAS',x[key+'_source'],title,'whole blood',x[key+'_fusion_p'],'FUSION Z '+x[key+'_fusion_z']+'; FDR '+x[key+'_fusion_fdr']+'; source GWAS overlap likely; gene-model association not causal')
 edge('gene:ENSG00000138722','phenotype:published_insomnia','published_whole_blood_FUSION_TWAS','https://doi.org/10.1152/jn.00490.2024','insomnia','whole blood','5.48e-5','MMRN1, FUSION Z -4.03; FDR 0.01; gene-model association, not rs7691487-specific or causal')
cp=pathlib.Path('outputs/stage4_cd4_smad5_coloc_abf_sensitivity.csv')
if cp.exists():
 for x in csv.DictReader(cp.open()):
  if float(x['p12'])==1e-5 and float(x['qtl_prior_sd'])==.15:
   edge('locus:chr5:135702653-136702653_GRCh38','gene:ENSG00000113658','exploratory_cell_type_region_colocalization','eQTL Catalogue r7 QTD000612 and '+x['trait'],x['trait'],'naive CD4 T cells','','H1='+x['H1']+'; H3='+x['H3']+'; H4='+x['H4']+'; single signal model; incomplete molecular credible-set coverage, priors and sample overlap uncertain; no causal inference')
edge('gene:ENSG00000113658','phenotype:published_insomnia','published_whole_blood_SMR','https://pmc.ncbi.nlm.nih.gov/articles/PMC12651937/','sleeplessness / insomnia','whole blood','3.25e-11','SMAD5 blood eQTLGen SMR reported in 2025 paper; exact insomnia GWAS and HEIDI for this association not independently verified; possible sample overlap; not variant-specific proof, CD4-specific evidence, frailty association or causal mediation')
ip=pathlib.Path('outputs/stage3_inka1_muscle_coloc_sensitivity.csv')
if ip.exists():
 for x in csv.DictReader(ip.open()):
  if float(x['p12'])==1e-5 and float(x['prior_scale'])==1:
   edge('locus:chr3:48803261-50803261_GRCh38','gene:ENSG00000185614','exploratory_region_colocalization','eQTL Catalogue r7 QTD000281 and '+x['trait'],x['trait'],'skeletal muscle','','H1='+x['H1']+'; H3='+x['H3']+'; H4='+x['H4']+'; 2-Mb region; single signal; incomplete molecular credible-set coverage and uncertain phenotype scales; no causal inference')
sp=pathlib.Path('outputs/stage4_mouse_brain_insomnia_scdrs_groups.csv')
if sp.exists():
 for x in csv.DictReader(sp.open()):
  edge('cell_type:Zeisel_mouse_'+x['group'],'phenotype:published_insomnia','exploratory_cross_species_cell_group_scDRS','Zeisel 2015 mouse brain and scDRS Jansen2019 insomnia gene set','insomnia','mouse cortex / hippocampal CA1',x['assoc_mcp'],'group association BH q across seven groups '+x['assoc_bh_fdr_7_groups']+'; N cells '+x['n_cell']+'; cross-species gene-program evidence; no frailty, direct variant, cell causality or donor-aware analysis')
rp=pathlib.Path('outputs/stage3_rs323509_testis_splice_credible_set.csv')
if rp.exists():
 for x in csv.DictReader(rp.open()):
  if x['rsid']=='rs323509':
   edge('variant:rs323509','gene:ENSG00000251574','bulk_splicing_QTL_credible_set','eQTL Catalogue r7 QTD000340',x['molecular_trait_id'],'testis',x['pvalue'],'molecular PIP '+x['pip']+'; set '+x['cs_id']+'; exact junction; related GTEx releases not independent; no disease mediation established')
hp=pathlib.Path('outputs/stage4_human_musc_insomnia_scdrs_groups.csv')
if hp.exists():
 for x in csv.DictReader(hp.open()):
  edge('cell_type:Kedlian_human_'+x['group'],'phenotype:published_insomnia','exploratory_human_cell_group_scDRS','Kedlian human muscle stem-cell subset and scDRS Jansen2019 insomnia gene set','insomnia','human skeletal muscle stem cells',x['assoc_mcp'],'group association BH q across four subtypes '+x['assoc_bh_fdr_4_subtypes']+'; N cells '+x['n_cell']+'; 12 represented donors; group test not donor-aware; no frailty, direct variant, or causal cell mechanism')
ep=pathlib.Path('outputs/stage4_human_musc_candidate_gene_expression.csv')
if ep.exists():
 for x in csv.DictReader(ep.open()):
  if int(x['donors_with_detected_expression'])>0:
   edge('gene:'+x['gene_id'],'cell_type:Kedlian_human_'+x['cell_group'],'descriptive_cell_type_gene_expression','Kedlian human muscle stem-cell subset',x['gene_symbols'],'human skeletal muscle stem cells','','equal-donor mean normalized expression '+x['equal_donor_mean_normalized_expression']+'; median donor fraction detected '+x['median_donor_fraction_detected']+'; detected in '+x['donors_with_detected_expression']+'/'+x['n_donors_at_least_10_cells']+' donors with >=10 cells; expression only, no genetic association or mechanism')
for path,locus,gene,tissue,source in [('outputs/stage3_rs323509_testis_coloc_sensitivity.csv','chr5:104246478-105246478_GRCh38','ENSG00000251574','testis','QTD000340 splicing'),('outputs/stage3_foxp2_pituitary_coloc_sensitivity.csv','chr7:113086327-115086327_GRCh38','ENSG00000128573','pituitary','QTD000301 expression')]:
 if pathlib.Path(path).exists():
  for x in csv.DictReader(open(path)):
   if float(x['p12'])==1e-5 and float(x['prior_scale'])==1:
    edge('locus:'+locus,'gene:'+gene,'exploratory_region_colocalization','eQTL Catalogue r7 '+source+' and '+x['trait'],x['trait'],tissue,'','H1='+x['H1']+'; H3='+x['H3']+'; H4='+x['H4']+'; single-signal model, prior-sensitive, phenotype scales uncertain; '+('high-PIP indel omitted' if gene=='ENSG00000128573' else 'testis relevance unknown')+'; no causal inference')
ap=pathlib.Path('outputs/stage3_adarb1_muscle_splice_coloc_sensitivity.csv')
if ap.exists():
 for x in csv.DictReader(ap.open()):
  if float(x['p12'])==1e-5 and float(x['prior_scale'])==1:
   edge('locus:chr21:44622459-45622459_GRCh38','gene:ENSG00000197381','exploratory_region_splice_colocalization','eQTL Catalogue r7 QTD000285 and '+x['trait'],x['trait'],'skeletal muscle','','H1='+x['H1']+'; H3='+x['H3']+'; H4='+x['H4']+'; single-signal, prior sensitive; 13/15 molecular credible variants; no causal inference')
op=pathlib.Path('outputs/stage3_original_gwas_coloc_sensitivity.csv')
if op.exists():
 info={'rs323509_testis':('chr5:104246478-105246478_GRCh38','ENSG00000251574','testis','QTD000340 splice'),'foxp2_pituitary':('chr7:113086327-115086327_GRCh38','ENSG00000128573','pituitary','QTD000301 expression'),'inka1':('chr3:48803261-50803261_GRCh38','ENSG00000185614','skeletal muscle','QTD000281 expression'),'adarb1_muscle_splice':('chr21:44622459-45622459_GRCh38','ENSG00000197381','skeletal muscle','QTD000285 splice')}
 for x in csv.DictReader(op.open()):
  if float(x['p12'])==1e-5 and float(x['prior_scale'])==1:
   locus,gene,tissue,qtl=info[x['analysis']]
   edge('locus:'+locus,'gene:'+gene,'exploratory_region_colocalization','eQTL Catalogue r7 '+qtl+' and '+x['trait'],x['trait'],tissue,'','H1='+x['H1']+'; H3='+x['H3']+'; H4='+x['H4']+'; pairwise single-signal, prior sensitive, candidate-selected; phenotype overlap and missing variants; no trait mediation established')
bp=pathlib.Path('outputs/stage3_remaining_coloc_baseline_summary.csv')
if bp.exists():
 import json
 ex={x['candidate']:x for x in json.load(open('outputs/stage3_remaining_catalogue_extract_receipt.json'))['items']}
 work={x['candidate']:x for x in csv.DictReader(open('outputs/stage3_remaining_molecular_worklist.csv')) if x['candidate']}
 for c,x in work.items():
  if ex[c].get('candidate_rows',0)==1:
   qpath=pathlib.Path('work/stage3_remaining_qtl')/(c+'.tsv');pos=x['marker_grch38'].split(':')[1]
   z=[v for v in csv.DictReader(qpath.open(),delimiter='\t') if v['rsid']==c and v['position']==pos and v['ref']==x['candidate_ref'] and v['alt']==x['candidate_alt']];assert len(z)==1
   edge('variant:'+c,'gene:'+x['gene_id'],'bulk_expression_QTL_reprocessing','eQTL Catalogue r7 '+x['dataset_r7'],'expression',x['tissue_v10'],z[0]['pvalue'],'nominal exact-marker QTL; related GTEx v8/v10 sources, not independent replication; no causality')
 for x in csv.DictReader(bp.open()):
  e=ex[x['candidate']];loc='chr'+e['chromosome']+':'+str(e['region'][0])+'-'+str(e['region'][1])+'_GRCh38'
  edge('locus:'+loc,'gene:'+x['gene_id'],'exploratory_region_colocalization','eQTL Catalogue r7 '+x['dataset']+' and '+x['trait'],x['trait'],x['tissue'],'','H0='+x['H0']+'; H1='+x['H1']+'; H2='+x['H2']+'; H3='+x['H3']+'; H4='+x['H4']+'; '+x['limitation']+'; credible set '+x['credible_set_coverage']+'; postselected, prior-sensitive, not causal')
np=pathlib.Path('outputs/stage3_negr1_coloc_sensitivity.csv')
if np.exists():
 edge('variant:rs2815764','gene:ENSG00000172260','bulk_expression_QTL_reprocessing','eQTL Catalogue r7 QTD000156','expression','brain caudate','5.63015e-05','nominal exact marker eQTL; v8/v10 GTEx share source, not independent replication; no causality')
 for x in csv.DictReader(np.open()):
  if float(x['p12'])==1e-5 and float(x['prior_scale'])==1:
   edge('locus:chr1:71288440-73288440_GRCh38','gene:ENSG00000172260','exploratory_region_colocalization','eQTL Catalogue r7 QTD000156 and '+x['trait'],x['trait'],'brain caudate','','H0='+x['H0']+'; H1='+x['H1']+'; H2='+x['H2']+'; H3='+x['H3']+'; H4='+x['H4']+'; fallback selection; no reported molecular credible set, prior-sensitive single-signal model; no causal inference')
mp=pathlib.Path('outputs/stage4_magma_cell_type_exploratory.csv')
if mp.exists():
 for x in csv.DictReader(mp.open()):
  edge('cell_type:FUMA:'+x['dataset']+'::'+x['cell_type'],'phenotype:'+x['trait'],'MAGMA_cell_type_gene_property','MAGMA v1.10 '+x['dataset'],x['trait'],x['domain'],x['p'],'BH q across six matrices for trait/model '+x['bh_q_all_six_matrices']+'; beta '+x['beta']+'; '+x['species']+'; '+('Brown 1-df gene model' if 'brown' in x['trait'] else 'SNPwise-mean gene model')+'; conditioned on Average expression; '+('regression failed, untestable; collinear with conditioned covariates; ' if not x['p'] else '')+'postselected matrix family; no specific gene or cell causality')
gp=pathlib.Path('outputs/stage4_ye_ffs_magma_gene_results.csv')
if gp.exists():
 ids={x['gene ID'].split('.')[0] for x in csv.DictReader(open('outputs/stage3_initial_variant_gene_table.csv')) if x['gene ID'].startswith('ENSG')}
 for x in csv.DictReader(gp.open()):
  if x['GENE'] in ids:
   edge('gene:'+x['GENE'],'phenotype:ye_ffs','MAGMA_gene_association','MAGMA v1.10 SNPwise-mean and Ye FFS','Ye FFS','',x['P'],'BH q across 19015 genes '+str(x['bh_q_all_19015_genes'])+'; Ensembl v92 symbol '+x['gene_symbol_v92']+'; LD-aware gene association, not a variant or mediation claim')
mfi=pathlib.Path('outputs/stage4_mak_fi_brown_magma_gene_results.csv')
if mfi.exists():
 ids={x['gene ID'].split('.')[0] for x in csv.DictReader(open('outputs/stage3_initial_variant_gene_table.csv')) if x['gene ID'].startswith('ENSG')}
 for x in csv.DictReader(mfi.open()):
  if x['GENE'] in ids:
   edge('gene:'+x['GENE'],'phenotype:mak_fi_brown','MAGMA_gene_association','MAGMA v1.10 Brown 1-df and Mak FI','Mak FI','',x['P'],'BH q across 19129 genes '+str(x['bh_q_all_genes'])+'; Ensembl v92 symbol '+x['gene_symbol_v92']+'; approximate LD-aware gene association, not a variant or mediation claim')
mf=pathlib.Path('outputs/stage4_mak_fi_magma_gene_results.csv')
if mf.exists():
 ids={x['gene ID'].split('.')[0] for x in csv.DictReader(open('outputs/stage3_initial_variant_gene_table.csv')) if x['gene ID'].startswith('ENSG')}
 for x in csv.DictReader(mf.open()):
  if x['GENE'] in ids:
   edge('gene:'+x['GENE'],'phenotype:mak_fi','MAGMA_gene_association','MAGMA v1.10 SNPwise-mean and Mak FI','Mak FI','',x['P'],'BH q across 19129 genes '+str(x['bh_q_all_genes'])+'; Ensembl v92 symbol '+x['gene_symbol_v92']+'; LD-aware gene association, not a variant or mediation claim')
for trait,label in [('ye_ffs_brown','Ye FFS'),('jansen_insomnia_brown','Jansen insomnia'),('atkins_fi_brown','Atkins FI')]:
 p=pathlib.Path('outputs/stage4_'+trait+'_magma_gene_results.csv')
 if p.exists():
  ids={x['gene ID'].split('.')[0] for x in csv.DictReader(open('outputs/stage3_initial_variant_gene_table.csv')) if x['gene ID'].startswith('ENSG')}
  for x in csv.DictReader(p.open()):
   if x['GENE'] in ids:
    edge('gene:'+x['GENE'],'phenotype:'+trait,'MAGMA_gene_association','MAGMA v1.10 Brown 1-df and '+label,label,'',x['P'],'BH q across all genes '+str(x['bh_q_all_genes'])+'; Ensembl v92 symbol '+x['gene_symbol_v92']+'; approximate LD-aware gene association, not a variant or mediation claim')
abc=pathlib.Path('outputs/stage3_abc_candidate_regulatory_links.csv')
if abc.exists():
 for x in csv.DictReader(abc.open()):
  edge('variant:'+x['candidate'],'gene_symbol:'+x['TargetGene'],'predicted_enhancer_gene_link','Nasser 2021 ABC thresholded predictions',x['TargetGene'],x['CellType'],'','ABC score '+x['ABC.Score']+'; '+x['class']+' interval '+x['chr']+':'+x['start']+'-'+x['end']+' GRCh37; exact marker position overlaps; model prediction only, no causal or allele-specific inference')
pf=pathlib.Path('outputs/stage3_positional_fallback_coloc_sensitivity.csv')
if pf.exists():
 info={'rs9845387':('ENSG00000185565','chr3:115707088-117707088_GRCh38','0.102125'),'rs4886140':('ENSG00000139734','chr13:58259385-60259385_GRCh38','0.580882')}
 for c,(gid,loc,prob) in info.items():
  edge('variant:'+c,'gene:'+gid,'nominal_QTL_evaluation','eQTL Catalogue r7 QTD000171','expression','brain cortex',prob,'positional fallback; exact-marker nominal p not significant; no credible set reported, no causal inference')
 for x in csv.DictReader(pf.open()):
  if float(x['p12'])==1e-5 and float(x['prior_scale'])==1:
   gid,loc,prob=info[x['candidate']]
   edge('locus:'+loc,'gene:'+gid,'exploratory_region_colocalization','eQTL Catalogue r7 QTD000171 and '+x['trait'],x['trait'],'brain cortex','','H0='+x['H0']+'; H1='+x['H1']+'; H2='+x['H2']+'; H3='+x['H3']+'; H4='+x['H4']+'; positional fallback; no molecular credible set, prior-sensitive single-signal model; no causal inference')
so=pathlib.Path('outputs/stage3_smim32_original_gwas_coloc_sensitivity.csv')
if so.exists():
 for x in csv.DictReader(so.open()):
  if float(x['p12'])==1e-5 and float(x['prior_scale'])==1:
   edge('locus:chr5:135702653-136702653_GRCh38','gene:ENSG00000271824','exploratory_region_colocalization','eQTL Catalogue r7 QTD000281 and '+x['trait'],x['trait'],'skeletal muscle','','H1='+x['H1']+'; H3='+x['H3']+'; H4='+x['H4']+'; only 74/98 molecular credible variants; overlapping samples and single signal; prior-sensitive, not causal')
fp=pathlib.Path('outputs/stage4_fuma_candidate_gene_expression.csv')
if fp.exists():
 for x in csv.DictReader(fp.open()):
  if float(x['processed_expression'])>0:
   edge('gene:'+x['gene_id'],'cell_type:FUMA:'+x['dataset']+'::'+x['cell_type'],'descriptive_cell_type_gene_expression','FUMA scRNA processed '+x['dataset'],x['gene_symbols'],x['domain'],'','source processed expression '+x['processed_expression']+'; '+x['species']+'; average expression, not a genetic association or causal link; no donor-aware inference')
tf=pathlib.Path('outputs/stage3_published_fi_twas_supplement_genes.csv')
if tf.exists():
 import collections
 ids_by_symbol=collections.defaultdict(set)
 for z in csv.DictReader(open('outputs/stage3_initial_variant_gene_table.csv')):
  if z['gene ID'].startswith('ENSG'):ids_by_symbol[z['Candidate gene']].add(z['gene ID'].split('.')[0])
 for x in csv.DictReader(tf.open()):
  if not x['project_candidate_markers']:continue
  ids=ids_by_symbol[x['gene_symbol']];node='gene:'+next(iter(ids)) if len(ids)==1 else 'gene_symbol:'+x['gene_symbol']
  edge(node,'phenotype:published_FI','published_FI_supplement_TWAS','Lin et al. 2024 supplementary table','FI','GTEx whole blood' if x['method'].startswith('FUSION') else 'GTEx cross-tissue',x['p'],'reported FDR '+x['reported_FDR']+'; gene symbol match; not marker-specific, independent replication or causal mediation')
pq=pathlib.Path('outputs/stage3_selected_molecular_qtl_ld_proxies.csv')
if pq.exists():
 for x in csv.DictReader(pq.open()):
  proxy='variant_GRCh37:'+x['proxy_GRCh37_chr']+':'+x['proxy_GRCh37_pos']+':'+x['proxy_allele_pair']
  edge(proxy,'gene:'+x['gene_id'],'selected_proxy_molecular_QTL','eQTL Catalogue r7 '+x['dataset'],x['molecular_trait_id'],x['tissue'],x['qtl_p_unadjusted'],'reference proxy '+x['proxy_reference_id']+' for '+x['candidate']+' r² '+x['candidate_proxy_r2']+'; source GRCh38 variant '+x['qtl_variant_grch38']+'; nominal selected-gene QTL; rsID/allele match across assumed builds; no independent liftover, multiple-testing correction or causal inference')
hm=pathlib.Path('outputs/stage4_hpa_muscle_candidate_expression.csv')
if hm.exists():
 for x in csv.DictReader(hm.open()):
  if float(x['nCPM'])>0:
   edge('gene:'+x['gene_id'],'cell_type:HPA_GSE233882:'+x['cluster'],'descriptive_cell_type_gene_expression','Human Protein Atlas GSE233882 processed cluster',x['candidate_gene_symbols'],'human skeletal muscle','','HPA cell type '+x['cell_type']+'; detail '+x['cell_type_detail']+'; '+x['cell_count']+' nuclei; nCPM '+x['nCPM']+'; source study included peripheral artery disease and controls; expression only; no genetic or causal inference')
hg=pathlib.Path('outputs/stage4_hpa_muscle_four_trait_cell_properties.csv')
if hg.exists():
 meta={x['Cluster']:x for x in csv.DictReader(open('outputs/stage4_hpa_skeletal_muscle_cluster_metadata.csv'))}
 for x in csv.DictReader(hg.open()):
  m=meta[x['cell_type_source_label']]
  edge('cell_type:HPA_GSE233882:'+x['cell_type_source_label'],'phenotype:'+x['trait'],'MAGMA_secondary_cell_type_gene_property','MAGMA v1.10 HPA GSE233882',x['trait'],'human skeletal muscle',x['p'],'BH q across 100 planned trait-cluster tests '+x['bh_q_100_trait_cell_tests']+'; source label '+m['Cell type detail']+'; '+m['Cell count']+' nuclei; source study included peripheral artery disease and controls; beta '+x['beta']+'; Brown 1-df; conditioned on equal-cluster mean log1p(nCPM); post hoc; no cell or gene causality'+('; untestable regression' if not x['p'] else ''))
pfc=pathlib.Path('outputs/stage4_adult_pfc_four_trait_cell_properties.csv')
if pfc.exists():
 for x in csv.DictReader(pfc.open()):
  edge('cell_type:FUMA:'+x['dataset']+'::'+x['cell_type_source_label'],'phenotype:'+x['trait'],'MAGMA_secondary_cell_type_gene_property','MAGMA v1.10 '+x['dataset'],x['trait'],'human prefrontal cortex',x['p'],'BH q across 60 planned trait-cell tests '+x['bh_q_60_trait_cell_tests']+'; beta '+x['beta']+'; Brown 1-df gene model; conditioned on Average; post hoc matrix; no specific gene or cell causality')
ldp=pathlib.Path('outputs/stage3_candidate_reference_ld_proxies.csv')
if ldp.exists():
 for x in csv.DictReader(ldp.open()):
  if x['proxy_pos_GRCh37']==x['candidate_pos_GRCh37']:continue
  proxy='variant_GRCh37:'+x['candidate_chr_GRCh37']+':'+x['proxy_pos_GRCh37']+':'+x['proxy_bim_A1_A2']
  edge('variant:'+x['candidate'],proxy,'reference_LD_proxy','1000 Genomes EUR MAGMA reference','LD','','','proxy rsID '+x['proxy_reference_id']+'; r² '+x['r2_genotype_dosage']+'; n='+x['n_pairwise_nonmissing']+'; 250 kb window; reference dosage correlation, not causality or GWAS fine mapping')
abcp=pathlib.Path('outputs/stage3_abc_ld_proxy_regulatory_links.csv')
if abcp.exists():
 mapping={(x['candidate'],x['abc_target_symbol']):x['candidate_gene_id_v92'] for x in csv.DictReader(open('outputs/stage3_abc_ld_proxy_gene_id_mapping.csv'))}
 for x in csv.DictReader(abcp.open()):
  proxy='variant_GRCh37:'+x['proxy_chr_GRCh37']+':'+x['proxy_pos_GRCh37']+':'+x['proxy_bim_A1_A2']
  gid=mapping[(x['candidate'],x['TargetGene'])];target='gene:'+gid if gid else 'gene_symbol:'+x['TargetGene']
  edge(proxy,target,'predicted_proxy_regulatory_link','Nasser 2021 ABC thresholded predictions',x['TargetGene'],x['CellType'],'','proxy '+x['proxy_reference_id']+' linked to candidate '+x['candidate']+' with reference r² '+x['proxy_r2']+'; '+x['class']+' interval '+x['chr']+':'+x['start']+'-'+x['end']+'; ABC score '+x['ABC.Score']+'; target '+('matched to Ensembl v92 by symbol, chromosome and TSS among tested coding genes' if gid else 'symbol unresolved')+'; model prediction only, no allele-specific or causal inference')
p=pathlib.Path('outputs/stage4_variant_gene_phenotype_network_scaffold.csv')
with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
print('edges',len(rows),'GWAS',sum(x['edge_type']=='GWAS_marginal_association' for x in rows),'molecular',sum(x['edge_type'].startswith('bulk_tissue_') for x in rows),'cell-type edges',sum('cell_type' in x['edge_type'] for x in rows))

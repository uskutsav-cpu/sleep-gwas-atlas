import csv,collections,math,pathlib
m={r['variant']:r for r in csv.DictReader(open('outputs/shared22_eight_source_evidence_matrix.csv'))}; extra=next(x for x in csv.DictReader(open('outputs/integrated_15_candidate_review.csv')) if x['variant']=='rs323509')
m['rs323509']={'variant':'rs323509','chr_GRCh38':extra['GRCh37_chr'],'position_GRCh38':extra['GRCh38_pos'],'score_counted_allele':extra['counted_allele'],'other_allele':extra['other_allele'],'nearest_gene_body_positional':extra['nearest_gene_body_positional'],'atkins_FI_beta_aligned':extra['FI_beta_on_counted_allele'],'atkins_FI_p':extra['FI_p'],'ye_FFS_beta_aligned':extra['FFS_beta_inferred_on_counted_allele'],'ye_FFS_p':extra['FFS_p'],'ukb_insomnia_beta_aligned':extra['UKB_insomnia_beta_on_counted_allele'],'ukb_insomnia_p':extra['UKB_insomnia_p']}
q=list(csv.DictReader(open('outputs/stage3_gtex_v10_candidate_qtls.csv'))); regpath=__import__('pathlib').Path('outputs/stage3_ensembl_regulatory_10kb.csv'); reg=list(csv.DictReader(regpath.open())) if regpath.exists() else [];groups=collections.defaultdict(list)
for r in q: groups[(r['candidate'],r['gene_symbol'],r['gencode_id'],r['tissue'])].append(r)
result=[]
for variant,x in m.items():
 keys=[key for key in groups if key[0]==variant]
 if not keys:keys=[(variant,x['nearest_gene_body_positional'],'','not returned in GTEx v10 QTL query')]
 for key in keys:
  rows=groups.get(key,[]); eq=[r for r in rows if r['qtl_type']=='Eqtl'];sq=[r for r in rows if r['qtl_type']=='Sqtl']
  def evidence(rs):
   if not rs:return 'No association returned for this exact gene/tissue/query'
   return 'GTEx v10 precomputed significant association; '+str(len(rs))+' row(s); min source p='+format(min(float(r['p_value']) for r in rs),'.4g')
  row={'Variant':variant,'Locus':f"chr{x['chr_GRCh38']}:{x['position_GRCh38']} (GRCh38 assumed mapping)",'Insomnia effect':f"UKB beta {x['ukb_insomnia_beta_aligned']} p={x['ukb_insomnia_p']}",'FI effect':f"Atkins beta {x['atkins_FI_beta_aligned']} p={x['atkins_FI_p']}",'FFS effect':f"Ye beta {x['ye_FFS_beta_aligned']} p={x['ye_FFS_p']}",'Candidate gene':key[1],'Tissue':key[3],'eQTL':evidence(eq),'sQTL':evidence(sq),'Colocalization':'Not tested','Evidence status':'Provisional exact-marker molecular association; shared causal signal not established; causality, independent signal and allele semantics unverified' if rows else 'No association returned; nearest positional gene only; not a negative result','counted allele / other':x['score_counted_allele']+'/'+x['other_allele'],'gene ID':key[2]}
  if variant=='rs7734045' and key[1]=='SMIM32' and key[3]=='Muscle_Skeletal':
   cp=__import__('pathlib').Path('outputs/stage3_smim32_finngen_coloc_abf_sensitivity.csv')
   if cp.exists():
    sensitivity=list(csv.DictReader(cp.open())); selected=[z for z in sensitivity if float(z['qtl_prior_effect_sd'])==.15 and float(z['p12'])==1e-5]; assert len(selected)==1
    row['Colocalization']='Exploratory region-level FinnGen insomnia vs SMIM32 eQTL: H4 '+format(float(selected[0]['H4']),'.4f')+', H1 '+format(float(selected[0]['H1']),'.3f')+' at illustrative priors; see report'
    row['Evidence status']+='; credible-set coverage incomplete, so single-signal coloc cannot support shared/distinct inference; not a marker-level causal claim'
   mp=__import__('pathlib').Path('outputs/stage3_smim32_mak_coloc_abf_sensitivity.csv')
   if mp.exists():
    mm=list(csv.DictReader(mp.open()));
    for trait in ('Mak_FI','Mak_FP'):
     selected=[z for z in mm if z['trait']==trait and float(z['gwas_prior_scale_factor'])==1 and float(z['p12'])==1e-5];assert len(selected)==1
     row['Colocalization']+='; '+trait+' H4 '+format(float(selected[0]['H4']),'.4f')+' H3 '+format(float(selected[0]['H3']),'.4f')
  if variant=='rs7734045' and key[1]=='SMAD5' and key[3]=='Muscle_Skeletal':
   cp=__import__('pathlib').Path('outputs/stage3_smad5_splice_coloc_abf_sensitivity.csv')
   if cp.exists():
    data=list(csv.DictReader(cp.open())); selected=[z for z in data if float(z['qtl_prior_effect_sd'])==.15 and float(z['p12'])==1e-5]; assert len(selected)==3
    row['Colocalization']='Exploratory region-level SMAD5 splice: '+ '; '.join(z['trait']+' H4 '+format(float(z['H4']),'.4f')+' H3 '+format(float(z['H3']),'.4f') for z in selected)+'; see report'
    row['Evidence status']+='; two molecular signals and incomplete credible-set coverage invalidate definitive single-signal inference; not a marker-level causal claim'
  row['Molecular fine mapping'] = 'not yet evaluated'
  if variant=='rs2838788' and key[3]=='Muscle_Skeletal' and key[1]=='ADARB1' and pathlib.Path('outputs/stage3_adarb1_muscle_splice_finemapping.md').exists():
   row['Molecular fine mapping']='Catalogue r7 muscle ADARB1 selected splice junction: one 15-variant credible set; candidate not listed; see report'
   if pathlib.Path('outputs/stage3_adarb1_muscle_splice_coloc_sensitivity.csv').exists():
    z=[a for a in csv.DictReader(open('outputs/stage3_adarb1_muscle_splice_coloc_sensitivity.csv')) if float(a['prior_scale'])==1 and float(a['p12'])==1e-5];assert len(z)==4
    row['Colocalization']='Exploratory ADARB1 splice region: '+ '; '.join(a['trait']+' H1 '+format(float(a['H1']),'.3f')+' H3 '+format(float(a['H3']),'.3f')+' H4 '+format(float(a['H4']),'.3f') for a in z)+'; prior sensitive and incomplete coverage; see report'
    row['Evidence status']+='; single-signal model and missing credible-set variants limit inference'
  if variant=='rs323509' and key[3]=='Testis' and key[2].split('.')[0]=='ENSG00000251574' and pathlib.Path('outputs/stage3_rs323509_testis_finemapping.md').exists():
   row['Molecular fine mapping']='Catalogue r7 testis junction chr5:104739412-104739497: 4-variant credible set; candidate molecular PIP 0.291; see report'
   if pathlib.Path('outputs/stage3_rs323509_testis_coloc_sensitivity.csv').exists():
    z=[a for a in csv.DictReader(open('outputs/stage3_rs323509_testis_coloc_sensitivity.csv')) if float(a['prior_scale'])==1 and float(a['p12'])==1e-5];assert len(z)==4
    row['Colocalization']='Exploratory testis splice region: '+ '; '.join(a['trait']+' H3 '+format(float(a['H3']),'.3f')+' H4 '+format(float(a['H4']),'.3f') for a in z)+'; prior-sensitive; see report'
    row['Evidence status']+='; tissue relevance unknown and single-signal model is prior-sensitive'
  if variant=='rs34290943' and key[3]=='Muscle_Skeletal' and key[1]=='INKA1' and pathlib.Path('outputs/stage3_inka1_muscle_finemapping.md').exists():
   row['Molecular fine mapping']='Catalogue r7 skeletal muscle INKA1: 1 reported 88-variant credible set; candidate not listed; see stage3_inka1_muscle_finemapping.md'
   if pathlib.Path('outputs/stage3_inka1_muscle_coloc_sensitivity.csv').exists():
    z=list(csv.DictReader(open('outputs/stage3_inka1_muscle_coloc_sensitivity.csv')));z=[a for a in z if float(a['prior_scale'])==1 and float(a['p12'])==1e-5];assert len(z)==4
    row['Colocalization']='Exploratory INKA1 muscle region: '+ '; '.join(a['trait']+' H1 '+format(float(a['H1']),'.3f')+' H3 '+format(float(a['H3']),'.3f')+' H4 '+format(float(a['H4']),'.4f') for a in z)+'; incomplete coverage; see report'
    row['Evidence status']+='; single-signal model and incomplete coverage limit inference'
  if key[3]=='Muscle_Skeletal' and key[1]=='SMIM32':
   row['Molecular fine mapping']='Catalogue r7 muscle expression: 1 credible set, 98 variants; candidate not listed; incomplete overlap with filtered GWAS inputs'
  if key[3]=='Muscle_Skeletal' and key[1]=='SMAD5':
   row['Molecular fine mapping']='Catalogue r7 muscle selected splice event: 2 credible sets, 52 and 18 variants; candidate not listed; incomplete GWAS coverage'
   if variant=='rs7734045':row['Colocalization']+='; FinnGen MAF 0.5% sensitivity H4 0.0890, L2 coverage 14/18; multiple signals remain'
  if key[3]=='Pituitary' and key[1]=='FOXP2' and variant in ('rs12536335','rs12671330','rs17312861'):
   row['Molecular fine mapping']='Catalogue r7 pituitary expression: 1 credible set, 4 variants; candidate not listed; indel missing from GWAS comparison'
   if pathlib.Path('outputs/stage3_foxp2_pituitary_coloc_sensitivity.csv').exists():
    z=[a for a in csv.DictReader(open('outputs/stage3_foxp2_pituitary_coloc_sensitivity.csv')) if float(a['prior_scale'])==1 and float(a['p12'])==1e-5];assert len(z)==4
    row['Colocalization']='Exploratory FOXP2 pituitary region: '+ '; '.join(a['trait']+' H3 '+format(float(a['H3']),'.3f')+' H4 '+format(float(a['H4']),'.3f') for a in z)+'; missing high-PIP indel; see report'
    row['Evidence status']+='; missing credible-set indel and single-signal model limit inference'
  sp=pathlib.Path('outputs/stage3_smim32_original_gwas_coloc_sensitivity.csv')
  if variant=='rs7734045' and key[1]=='SMIM32' and key[3]=='Muscle_Skeletal' and sp.exists():
   z=[a for a in csv.DictReader(sp.open()) if float(a['prior_scale'])==1 and float(a['p12'])==1e-5];assert len(z)==2
   row['Colocalization']+='; original GWAS '+ '; '.join(a['trait']+' H1 '+format(float(a['H1']),'.3f')+' H3 '+format(float(a['H3']),'.3f')+' H4 '+format(float(a['H4']),'.3f') for a in z)+'; 74/98 molecular variants matched; see report'
  np=pathlib.Path('outputs/stage3_negr1_coloc_sensitivity.csv')
  if variant=='rs2815764' and key[1]=='NEGR1' and key[3]=='Brain_Caudate_basal_ganglia' and np.exists():
   z=[a for a in csv.DictReader(np.open()) if float(a['prior_scale'])==1 and float(a['p12'])==1e-5];assert len(z)==6
   row['Molecular fine mapping']='Catalogue r7 caudate NEGR1: no reported molecular credible set; not a negative QTL result'
   row['Colocalization']='Exploratory regional model: '+ '; '.join(a['trait']+' H1 '+format(float(a['H1']),'.3f')+' H3 '+format(float(a['H3']),'.3f')+' H4 '+format(float(a['H4']),'.3f') for a in z if a['trait'] in ('Atkins_FI','Jansen_insomnia'))+'; other traits and H0/H2 in report; no reported credible set'
   row['Evidence status']+='; fallback selected after another gene had no v8 data; prior-sensitive and not causal'
  batchp=pathlib.Path('outputs/stage3_remaining_coloc_baseline_summary.csv')
  if batchp.exists():
   work=[a for a in csv.DictReader(open('outputs/stage3_remaining_molecular_worklist.csv')) if a['candidate']==variant and a['gene_id']==key[2].split('.')[0] and a['tissue_v10']==key[3]]
   if work:
    fin=next(a for a in csv.DictReader(open('outputs/stage3_remaining_molecular_finemapping_summary.csv')) if a['candidate']==variant)
    row['Molecular fine mapping']='Catalogue r7 '+fin['dataset']+': '+fin['credible_set_count']+' reported molecular set(s); '+fin['credible_set_sizes']+'; candidate membership '+fin['candidate_in_any_set']+'; see batch review'
    z=[a for a in csv.DictReader(batchp.open()) if a['candidate']==variant and a['trait'] in ('Atkins_FI','Jansen_insomnia')]
    if z:
     assert len(z)==2
     row['Colocalization']='Exploratory regional model: '+ '; '.join(a['trait']+' H1 '+format(float(a['H1']),'.3f')+' H3 '+format(float(a['H3']),'.3f')+' H4 '+format(float(a['H4']),'.3f') for a in z)+'; H0/H2 and other GWAS in batch review; '+z[0]['limitation']
     row['Evidence status']+='; postselected single-signal screening, prior-sensitive and not causal'
    elif fin['region_gene_rows']=='0':row['Evidence status']+='; no selected gene rows returned in related GTEx v8 Catalogue release'
  originalp=pathlib.Path('outputs/stage3_original_gwas_coloc_sensitivity.csv')
  targets={('rs323509','ENSG00000251574','Testis'):'rs323509_testis',('rs34290943','INKA1','Muscle_Skeletal'):'inka1',('rs2838788','ADARB1','Muscle_Skeletal'):'adarb1_muscle_splice'}
  analysis=targets.get((variant,key[1],key[3]))
  if key[1]=='FOXP2' and key[3]=='Pituitary' and variant in ('rs12536335','rs12671330','rs17312861'):analysis='foxp2_pituitary'
  if originalp.exists() and analysis:
   z=[a for a in csv.DictReader(originalp.open()) if a['analysis']==analysis and float(a['prior_scale'])==1 and float(a['p12'])==1e-5];assert len(z)==2
   text='original GWAS '+ '; '.join(a['trait']+' H3 '+format(float(a['H3']),'.3f')+' H4 '+format(float(a['H4']),'.3f') for a in z)+'; pairwise, see original GWAS report'
   row['Colocalization']=('Exploratory '+text) if row['Colocalization']=='Not tested' else (row['Colocalization']+'; '+text)
  nearby=[z for z in reg if z['variant']==variant]
  row['Ensembl regulatory features within 10kb'] = ('; '.join(z['regulatory_feature_id']+' '+z['description']+' '+z['distance_bp']+'bp' for z in sorted(nearby,key=lambda z:int(z['distance_bp']))) if nearby else ('none returned' if regpath.exists() else 'not queried'))
  row['regulatory_link_status'] = 'Positional only: no tissue-specific variant function or target gene established'
  result.append(row)
cellp=__import__('pathlib').Path('outputs/stage4_onek1k_candidate_credible_set_matches.csv')
if cellp.exists():
 seen=set()
 for x in csv.DictReader(cellp.open()):
  k=(x['candidate'],x['qtl_variant_id'],x['dataset'],x['molecular_trait_id'],x['credible_set_id'])
  if k in seen:continue
  seen.add(k); base=next(z for z in result if z['Variant']==x['candidate']).copy()
  base['Candidate gene']='SMAD5' if x['gene_id']=='ENSG00000113658' else x['gene_id'];base['gene ID']=x['gene_id'];base['Tissue']='CD4 naive T cells (OneK1K PBMC)'
  base['eQTL']='Catalogue r7 '+x['dataset']+' credible-set membership; PIP '+x['pip']+'; nominal source p '+x['p_value'];base['sQTL']='not tested in this cell-specific lookup';base['Colocalization']='Not tested for this cell-specific QTL'
  if x['gene_id']=='ENSG00000113658' and pathlib.Path('outputs/stage4_cd4_smad5_coloc_abf_sensitivity.csv').exists():
   cs=list(csv.DictReader(open('outputs/stage4_cd4_smad5_coloc_abf_sensitivity.csv')));cs=[z for z in cs if float(z['p12'])==1e-5 and float(z['qtl_prior_sd'])==.15];assert len(cs)==3
   base['Colocalization']='Exploratory CD4 region: '+ '; '.join(z['trait']+' H1 '+format(float(z['H1']),'.3f')+' H3 '+format(float(z['H3']),'.3f')+' H4 '+format(float(z['H4']),'.3f') for z in cs)+'; incomplete credible-set coverage; see report'
  base['Evidence status']='Pseudobulk cell-type molecular credible-set association; not disease colocalization, enrichment or proof of cell-specific causality'
  base['Molecular fine mapping']='Catalogue r7 '+x['dataset']+' credible set '+x['credible_set_id']+'; set size '+x['credible_set_size']+'; PIP '+x['pip']
  result.append(base)
if __import__('pathlib').Path('outputs/stage3_published_fi_insomnia_twas_review.md').exists():
 base=next(z for z in result if z['Variant']=='rs7691487' and z['Candidate gene']=='MMRN1').copy()
 base['Tissue']='Whole blood gene-prediction model (published insomnia FUSION)';base['eQTL']='Not directly tested for this variant and tissue';base['sQTL']='Not tested for this variant and tissue';base['Colocalization']='Not tested for this candidate versus published TWAS signal'
 base['Evidence status']='Published gene-level insomnia TWAS: FUSION z=-4.03, p=5.48e-5, FDR=0.01; the rs7691487–TWAS signal relationship is unverified and source GWAS may overlap'
 base['Molecular fine mapping']='Published FOCUS gene PIP has discrepant text (0.998) and table (0.849 for MMRN1); no resolution; not a variant PIP'
 result.append(base)
fallback=pathlib.Path('outputs/stage3_positional_fallback_coloc_sensitivity.csv')
if fallback.exists():
 for variant,gene,gid,prob in [('rs9845387','LSAMP','ENSG00000185565','0.102125'),('rs4886140','DIAPH3','ENSG00000139734','0.580882')]:
  base=next(z for z in result if z['Variant']==variant).copy();base['Candidate gene']=gene;base['gene ID']=gid;base['Tissue']='Brain_Cortex'
  base['eQTL']='Catalogue r7 cortex exact marker nominal p='+prob+'; not significant; positional fallback';base['sQTL']='Not tested for this fallback';base['Evidence status']='Positional gene selected without exact-marker v10 QTL; v8 cortex candidate p not significant; no causal target established'
  base['Molecular fine mapping']='No reported molecular credible set in Catalogue r7 cortex; absence is not negative evidence'
  z=[a for a in csv.DictReader(fallback.open()) if a['candidate']==variant and a['trait'] in ('Atkins_FI','Jansen_insomnia') and float(a['prior_scale'])==1 and float(a['p12'])==1e-5];assert len(z)==2
  base['Colocalization']='Exploratory positional fallback: '+ '; '.join(a['trait']+' H0 '+format(float(a['H0']),'.3f')+' H1 '+format(float(a['H1']),'.3f')+' H2 '+format(float(a['H2']),'.3f')+' H3 '+format(float(a['H3']),'.3f')+' H4 '+format(float(a['H4']),'.3f') for a in z)+'; no molecular credible set; see fallback report'
  result.append(base)
abc=pathlib.Path('outputs/stage3_abc_candidate_regulatory_links.csv')
if abc.exists():
 for a in csv.DictReader(abc.open()):
  base=next(z for z in result if z['Variant']==a['candidate']).copy();base['Candidate gene']=a['TargetGene'];base['gene ID']='';base['Tissue']='ABC prediction: '+a['CellType']
  base['eQTL']='Not established by ABC prediction';base['sQTL']='Not established by ABC prediction';base['Colocalization']='Not tested for this predicted target';base['Evidence status']='Predicted regulatory link from exact marker position to source target label; tissue and allele effect unproven'
  base['Molecular fine mapping']='Not evaluated for this predicted target';base['regulatory_link_status']='ABC '+a['class']+' interval '+a['chr']+':'+a['start']+'-'+a['end']+' GRCh37; score '+a['ABC.Score']+'; target symbol not independently resolved'
  result.append(base)
proxyabc=pathlib.Path('outputs/stage3_abc_ld_proxy_regulatory_links.csv')
if proxyabc.exists():
 from collections import defaultdict
 groups=defaultdict(list)
 for a in csv.DictReader(proxyabc.open()):groups[(a['candidate'],a['TargetGene'])].append(a)
 mapping={(x['candidate'],x['abc_target_symbol']):x for x in csv.DictReader(open('outputs/stage3_abc_ld_proxy_gene_id_mapping.csv'))}
 for (variant,target),group in sorted(groups.items()):
  base=next(z for z in result if z['Variant']==variant).copy();base['Candidate gene']=target;base['gene ID']=mapping[(variant,target)]['candidate_gene_id_v92']
  base['Tissue']='ABC predicted target across '+str(len({x['CellType'] for x in group}))+' source cell/tissue labels; see proxy link table'
  base['eQTL']='Not established by LD and ABC prediction';base['sQTL']='Not established by LD and ABC prediction';base['Colocalization']='Not tested for this predicted target'
  base['Evidence status']='Reference LD proxy overlaps predicted regulatory interval; '+mapping[(variant,target)]['mapping_status']+'; no allele effect or mechanism established'
  base['Molecular fine mapping']='LD proxies not necessarily GWAS credible variants; molecular fine mapping not evaluated for this predicted target'
  base['regulatory_link_status']='ABC '+','.join(sorted({x['class'] for x in group}))+'; '+str(len({(x['proxy_reference_id'],x['proxy_pos_GRCh37']) for x in group}))+' reference proxies; '+str(len(group))+' prediction records; '+mapping[(variant,target)]['mapping_status']
  result.append(base)
gfiles=[('Ye FFS MAGMA SNPwise-mean','outputs/stage4_ye_ffs_magma_gene_results.csv','bh_q_all_19015_genes'),('Mak FI MAGMA SNPwise-mean','outputs/stage4_mak_fi_magma_gene_results.csv','bh_q_all_genes'),('Mak FI MAGMA Brown-1df','outputs/stage4_mak_fi_brown_magma_gene_results.csv','bh_q_all_genes'),('Ye FFS MAGMA Brown-1df','outputs/stage4_ye_ffs_brown_magma_gene_results.csv','bh_q_all_genes'),('Jansen insomnia MAGMA Brown-1df','outputs/stage4_jansen_insomnia_brown_magma_gene_results.csv','bh_q_all_genes'),('Atkins FI MAGMA Brown-1df','outputs/stage4_atkins_fi_brown_magma_gene_results.csv','bh_q_all_genes')]
gmaps=[]
for label,path,qcol in gfiles:
 if pathlib.Path(path).exists():gmaps.append((label,{x['GENE']:x for x in csv.DictReader(open(path))},qcol))
for row in result:
 gid=row['gene ID'].split('.')[0];ev=[]
 for label,gm,qcol in gmaps:
  if gid in gm:ev.append(label+' gene p='+gm[gid]['P']+' BH q='+gm[gid][qcol])
 row['GWAS gene evidence']='; '.join(ev)+'; gene-level association only, not evidence for this variant or mechanism' if ev else 'No matching gene-based result available'
twaspath=pathlib.Path('outputs/stage3_published_fi_twas_supplement_genes.csv')
if twaspath.exists():
 twas=list(csv.DictReader(twaspath.open()))
 for row in result:
  hits=[x for x in twas if x['gene_symbol']==row['Candidate gene']]
  row['Published FI TWAS evidence']='; '.join(x['method']+' p='+x['p']+' reported FDR='+x['reported_FDR'] for x in hits)+'; gene-level symbol match, not marker-specific or causal' if hits else 'No match in source significant-gene tables; not a negative TWAS result'
assert len(m)==23 and len({r['Variant'] for r in result})==23
with open('outputs/stage3_initial_variant_gene_table.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=result[0].keys());w.writeheader();w.writerows(result)
print('23 marker coverage; gene/tissue rows',len(result),'placeholder-only',sum(not any(qr['candidate']==v for qr in q) for v in m))

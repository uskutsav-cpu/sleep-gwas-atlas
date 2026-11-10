#!/usr/bin/env python3
"""Build a traceable, non-training graph from the existing selected-candidate table."""
import csv,json,hashlib,math
from pathlib import Path
p=Path('outputs');src=p/'integrated_15_candidate_review.csv'
rows=list(csv.DictReader(src.open()));assert len(rows)==15 and len({r['variant'] for r in rows})==15
sha=hashlib.sha256(src.read_bytes()).hexdigest();src22=p/'expanded_22_candidate_risk_gene_table.csv'; rows22=list(csv.DictReader(src22.open())); assert len(rows22)==22 and len({r['variant'] for r in rows22})==22; sha22=hashlib.sha256(src22.read_bytes()).hexdigest(); ids22={r['variant'] for r in rows22}; nodes={};edges=[];variant_ids={}
def node(i,typ,**attrs):
 v={'id':i,'type':typ,**attrs}
 if i in nodes:assert nodes[i]==v
 nodes[i]=v
def number(x):
 try:
  y=float(x)
  return y if math.isfinite(y) else None
 except (ValueError,TypeError):return None
studies=[('FI_meta','FI_beta_on_counted_allele','FI_p','Atkins FI meta-analysis; includes UKB and TwinGene','FI effect','source effect direction cohort coding caveat'),('Ye_FFS','FFS_beta_inferred_on_counted_allele','FFS_p','Ye UK Biobank Fried frailty score','FFS effect','A1 effect convention inferred; cohort overlaps UKB'),('Jansen_UKB_insomnia','UKB_insomnia_beta_on_counted_allele','UKB_insomnia_p','Jansen UK Biobank insomnia symptoms','log odds','UKB overlap with frailty GWAS'),('FinnGen_insomnia','FinnGen_insomnia_beta_on_counted_allele','FinnGen_insomnia_p','FinnGen clinical insomnia','log odds','FinnGen study overlap possible'),('MVP_insomnia','MVP_insomnia_beta_on_counted_allele','MVP_insomnia_p','MVP clinical insomnia','log odds','sample overlap unverified'),('FinnGen_HFRS','FinnGen_HFRS_beta_on_counted_allele','FinnGen_HFRS_p','FinnGen hospital frailty risk score','beta','HFRS differs from FI/FFS; FinnGen overlap possible')]
for tag,beta,pval,desc,scale,caveat in studies:node('study:'+tag,'study_phenotype',description=desc,effect_scale=scale,caveat=caveat)
for r in rows:
 v='variant:GRCh37:'+r['GRCh37_chr']+':'+r['GRCh37_pos']+':'+r['counted_allele']+':'+r['other_allele'];node(v,'variant',rsid=r['variant'],build='GRCh37_inferred',chromosome=r['GRCh37_chr'],position=int(r['GRCh37_pos']),counted_allele=r['counted_allele'],other_allele=r['other_allele'],flags=r['flags'],candidate_lists=['15','22'] if r['variant'] in ids22 else ['15']);variant_ids[r['variant']]=v
 for name in r['nearest_gene_body_positional'].split(';'):
  if name:
   g='gene_symbol:'+name;node(g,'gene_symbol_unresolved',symbol=name);edges.append({'source':v,'target':g,'relation':'nearest_gene_body_positional','distance_bp':number(r['distance_bp']),'evidence':'positional_only_not_causal','input_file':src.name,'input_sha256':sha})
 for tag,beta,pval,desc,scale,caveat in studies:
  b=number(r[beta]);pv=number(r[pval]);
  if b is not None or pv is not None:
   assert pv is None or 0<=pv<=1
   edges.append({'source':v,'target':'study:'+tag,'relation':'gwas_marginal_association','effect_on_counted_allele':b,'p_value':pv,'input_file':src.name,'input_sha256':sha,'study_caveat':caveat,'postselected':True})
# Add 22-list-only variants in their source-table orientation; overlapping records use the existing 15-list orientation.
links22=[('FI_meta','FI_beta_aligned','FI_p'),('Ye_FFS','FFS_beta_aligned','FFS_p'),('Jansen_UKB_insomnia','UKB_insomnia_beta_aligned','UKB_insomnia_p'),('MVP_insomnia','MVP_insomnia_log_OR_aligned_from_rounded_OR','MVP_insomnia_p'),('FinnGen_HFRS','HFRS_beta_aligned','HFRS_p')]
for r in rows22:
 if r['variant'] in variant_ids:
  n=nodes[variant_ids[r['variant']]];assert n['chromosome']==r['GRCh37_chr_assumed'] and n['position']==int(r['GRCh37_pos_assumed']) and {n['counted_allele'],n['other_allele']}=={r['counted_allele'],r['other_allele']};continue
 v='variant:GRCh37:'+r['GRCh37_chr_assumed']+':'+r['GRCh37_pos_assumed']+':'+r['counted_allele']+':'+r['other_allele'];variant_ids[r['variant']]=v
 node(v,'variant',rsid=r['variant'],build='GRCh37_assumed',chromosome=r['GRCh37_chr_assumed'],position=int(r['GRCh37_pos_assumed']),counted_allele=r['counted_allele'],other_allele=r['other_allele'],flags='source_reference_alleles:'+r['reference_alleles'],candidate_lists=['22'])
 for name in r['nearest_gene_body_positional'].split(';'):
  if name:
   g='gene_symbol:'+name;node(g,'gene_symbol_unresolved',symbol=name);edges.append({'source':v,'target':g,'relation':'nearest_gene_body_positional','distance_bp':number(r['distance_bp']),'evidence':'positional_only_not_causal','input_file':src22.name,'input_sha256':sha22})
 for tag,beta,pval in links22:
  b=number(r[beta]);pv=number(r[pval]);assert pv is None or 0<=pv<=1
  if b is not None or pv is not None:edges.append({'source':v,'target':'study:'+tag,'relation':'gwas_marginal_association','effect_on_counted_allele':b,'p_value':pv,'input_file':src22.name,'input_sha256':sha22,'study_caveat':nodes['study:'+tag]['caveat']+'; MVP log OR derives from rounded source OR' if tag=='MVP_insomnia' else nodes['study:'+tag]['caveat'],'postselected':True})
assert len(variant_ids)==33
# Molecular credible sets are modeled molecular signals, not disease colocalizations.
for kind,file in [('eqtl','stage6_varshney_muscle_candidate_screen.json'),('caqtl','stage6_varshney_muscle_caqtl_candidate_screen.json')]:
 path=p/file;receipt=json.loads(path.read_text());receipt_sha=hashlib.sha256(path.read_bytes()).hexdigest()
 for m in receipt['matches']:
  if m['candidate'] not in variant_ids:continue
  assert m['pair_match'] is True and 0<=m['pip']<=1
  parts=m['molecular_set'].split('/');cell=parts[1];base=parts[-1];set_id='molecular_set:Varshney2026:'+m['molecular_set'];cell_id='cell:Varshney2026:'+cell
  node(set_id,'molecular_credible_set',assay=kind,archive_path=m['molecular_set'],cell_context=cell,tissue='human_skeletal_muscle',source=receipt['source'],source_archive_md5=receipt['credible_archive_md5'])
  node(cell_id,'cell_context',cell_label=cell,tissue='human_skeletal_muscle',source=receipt['source'])
  if not any(e['source']==set_id and e['target']==cell_id and e['relation']=='measured_in_cell_context' for e in edges):edges.append({'source':set_id,'target':cell_id,'relation':'measured_in_cell_context','input_file':file,'input_sha256':receipt_sha})
  edges.append({'source':variant_ids[m['candidate']],'target':set_id,'relation':'member_of_molecular_credible_set','molecular_pip':m['pip'],'variant_GRCh38_source':m['variant'],'allele_pair_match':m['pair_match'],'input_file':file,'input_sha256':receipt_sha,'caveat':'Molecular PIP only; not disease colocalization, causal gene, or association direction.'})
  if kind=='eqtl':
   gene=base.split('__L')[0].split('.')[-1];g='gene_symbol:'+gene;node(g,'gene_symbol_unresolved',symbol=gene)
   if not any(e['source']==set_id and e['target']==g and e['relation']=='expression_phenotype_gene' for e in edges):edges.append({'source':set_id,'target':g,'relation':'expression_phenotype_gene','gene_label_caveat':'expression phenotype, not established disease effector','input_file':file,'input_sha256':receipt_sha})
# One published immune-cell peak association and molecular–molecular coloc; not a disease-gene edge.
lookupfile=p/'stage6_tenk10k_cd4_peak_variant_lookup.json';lookup=json.loads(lookupfile.read_text());target=[r for r in lookup['results'] if r['variant_id']=='5:136202653:A:C' and r['phenotype_id']=='chr5:136132117-136134284'];assert len(target)==1 and variant_ids['rs7734045'] in nodes
rr=target[0];peak='peak:TenK10K:GRCh38:chr5:136132117-136134284';cell='cell:TenK10K:CD4_Naive';gene='gene_ensembl:ENSG00000113658'
node(peak,'chromatin_peak',build='GRCh38',interval='chr5:136132117-136134284',study='TenK10K',caveat='candidate variant lies outside this peak');node(cell,'cell_context',cell_label='CD4_Naive',tissue='peripheral_blood_mononuclear_cells',source='https://huggingface.co/datasets/anglixue/TenK10K_multiome');node(gene,'gene_ensembl',ensembl_id='ENSG00000113658',symbol='SMAD5')
sourcehash=hashlib.sha256(lookupfile.read_bytes()).hexdigest();edges.append({'source':variant_ids['rs7734045'],'target':peak,'relation':'caqtl_marginal_association','effect_allele_GRCh38':'C','variant_GRCh38':rr['variant_id'],'beta_per_alt':rr['slope'],'se':rr['slope_se'],'p_nominal':rr['pval_nominal'],'alt_frequency':rr['af'],'input_file':lookupfile.name,'input_sha256':sourcehash,'caveat':'Molecular association only; indexed source read, full source hash not verified; beta is not a disease effect.'});edges.append({'source':peak,'target':cell,'relation':'measured_in_cell_context','input_file':lookupfile.name,'input_sha256':sourcehash})
colfile=p/'stage6_tenk10k_smad5_molecular_coloc.csv';chosen=[r for r in csv.DictReader(colfile.open()) if r['peak']=='chr5:136132117-136134284' and r['celltype']=='CD4_Naive' and r['gene']=='ENSG00000113658'];assert len(chosen)==1
c=chosen[0];edges.append({'source':peak,'target':gene,'relation':'molecular_caqtl_eqtl_colocalization','PP_H4':number(c['PP.H4.abf']),'PP_H3':number(c['PP.H3.abf']),'tested_variants':int(c['nsnps_coloc_tested']),'source_top_variant_GRCh38':c['top_snp'],'input_file':colfile.name,'input_sha256':hashlib.sha256(colfile.read_bytes()).hexdigest(),'caveat':'Molecular–molecular, not disease colocalization. Source correction status provisional; model- and prior-dependent, donor overlap possible.'})
# Retain the model-based disease/molecular comparison as its own edge, including contrary evidence.
node('study:Mak_FI','study_phenotype',description='Mak UK Biobank frailty index',effect_scale='FI effect',caveat='Includes UKB; FI and insomnia definitions overlap; selected-locus analysis')
fitfile=p/'stage6_tenk10k_expanded_susie_coloc_sensitivity.csv'; fits=list(csv.DictReader(fitfile.open()));assert len(fits)==3
for f in fits:
 edges.append({'source':peak,'target':'study:Mak_FI','relation':'disease_molecular_colocalization_sensitivity','method':'SuSiE-coloc','nsnps':int(f['nsnps']),'p1':0.0001,'p2':0.0001,'p12':float(f['p12']),'PP_H3':float(f['PP.H3.abf']),'PP_H4':float(f['PP.H4.abf']),'qtl_signal_hit_GRCh38':f['hit1'],'fi_signal_hit_GRCh38':f['hit2'],'ld_reference':'1000 Genomes European 503 participants, rank-deficient','input_file':fitfile.name,'input_sha256':hashlib.sha256(fitfile.read_bytes()).hexdigest(),'caveat':'Postselected conditional model; LD mismatch, approximate sample sizes, reference mismatch and missing variants. Distinct-signal support is not proof of absence or a disease-gene effect.'})
# Selected FOXP2-region disease-disease association is a study pair, not a variant–gene causal edge.
foxfile=p/'stage6_foxp2_gnomad_ancestry_matched_coloc.csv';fox=[r for r in csv.DictReader(foxfile.open()) if r['mvp_n_assumption']=='effN'];assert len(fox)==3 and all(int(x['nsnps'])==1257 for x in fox)
region='region:FOXP2_selected_1Mb_GRCh37';node(region,'genomic_region',build='GRCh37',label='Selected one-megabase window around rs12536335',caveat='Postselected narrow window; does not include all regional signals')
for tag in ['MVP_insomnia','FinnGen_HFRS']:edges.append({'source':region,'target':'study:'+tag,'relation':'regional_gwas_support','input_file':foxfile.name,'input_sha256':hashlib.sha256(foxfile.read_bytes()).hexdigest(),'caveat':'Association signals evaluated in a postselected region; membership is not causal-gene evidence.'})
for x in fox:
 edges.append({'source':'study:MVP_insomnia','target':'study:FinnGen_HFRS','relation':'disease_disease_colocalization_sensitivity','region':region,'method':'SuSiE-coloc','nsnps':int(x['nsnps']),'p1':0.0001,'p2':0.0001,'p12':float(x['p12']),'PP_H3':float(x['PP.H3.abf']),'PP_H4':float(x['PP.H4.abf']),'insomnia_signal_hit':x['hit1'],'hfrs_signal_hit':x['hit2'],'ld_reference':'gnomAD NFE for MVP and FIN for FinnGen; external not in-sample','input_file':foxfile.name,'input_sha256':hashlib.sha256(foxfile.read_bytes()).hexdigest(),'caveat':'Postselected regional association, not gene or cell causality; sample overlap unverified, narrower window and external LD.'})
# Public aging-brain molecular associations; one mega-analysis per context avoids treating overlapping studies as replicates.
brainfile=p/'stage6_niagads_xqtl_candidate_screen.json';brain=json.loads(brainfile.read_text())
braindigest=hashlib.sha256(brainfile.read_bytes()).hexdigest()
brainrows=[x for x in brain['results']['rs7734045']['rows'] if x['xQTLtype']=='snuc-eQTL' and x['target']=='ENSG00000113658' and x['study']=='ROSMAP_CUIMC1_2_MIT' and x['variantId']=='chr5:136202653:A:C']
assert {r['context'] for r in brainrows}=={'exc','inh','oli'} and len(brainrows)==3
csfile=p/'stage6_niagads_smad5_cs_audit.json';cs=json.loads(csfile.read_text());csdigest=hashlib.sha256(csfile.read_bytes()).hexdigest()
for r in brainrows:
 ctx=r['context']; assay='assay:NIAGADS:ROSMAP_CUIMC1_2_MIT:SMAD5:'+ctx;cell='cell:NIAGADS:ROSMAP:'+ctx
 node(cell,'cell_context',cell_label={'exc':'excitatory neuron','inh':'inhibitory neuron','oli':'oligodendrocyte'}[ctx],tissue='human_brain',source='https://dss.niagads.org/datasets/NG00184/')
 node(assay,'molecular_association_context',gene='ENSG00000113658',assay='single_nucleus_eQTL',context=ctx,study='ROSMAP_CUIMC1_2_MIT',build='GRCh38',caveat='Mega-analysis includes overlapping ROSMAP sources; postselected; molecular association is not causal disease evidence.')
 edges.append({'source':variant_ids['rs7734045'],'target':assay,'relation':'single_nucleus_eqtl_marginal_association','variant_GRCh38':r['variantId'],'effect_allele_GRCh38':'C','beta_per_alt':float(r['beta']),'se':float(r['se']),'p_nominal':float(r['pvalue']),'input_file':brainfile.name,'input_sha256':braindigest,'caveat':'Atlas HMT-filtered; selected SNP-level molecular association; sign is not disease risk.'})
 edges.append({'source':assay,'target':cell,'relation':'measured_in_cell_context','input_file':brainfile.name,'input_sha256':braindigest})
 edges.append({'source':assay,'target':gene,'relation':'expression_phenotype_gene','input_file':brainfile.name,'input_sha256':braindigest,'gene_label_caveat':'Molecular phenotype, not disease effector.'})
 for match in cs['results'][ctx]['candidate']:
  assert match['variant_id']==r['variantId'];setid='molecular_set:NIAGADS:'+match['cs_95']
  node(setid,'molecular_credible_set',assay='single_nucleus_eQTL',cell_context=ctx,study='ROSMAP_CUIMC1_2_MIT',source='https://dss.niagads.org/datasets/NG00184/')
  edges.append({'source':variant_ids['rs7734045'],'target':setid,'relation':'member_of_molecular_credible_set','molecular_pip':match['pip'],'variant_GRCh38_source':match['variant_id'],'input_file':csfile.name,'input_sha256':csdigest,'caveat':'Conditional molecular PIP; not a disease colocalization probability.'})
  edges.append({'source':setid,'target':assay,'relation':'fine_maps_molecular_association_context','input_file':csfile.name,'input_sha256':csdigest})
assert len(cs['results']['oli']['candidate'])==0  # Missing from the selected set is not encoded as a negative edge.
brainfit=p/'stage6_niagads_smad5_susie_coloc_sensitivity.csv';brainfits=list(csv.DictReader(brainfit.open()));assert len(brainfits)==3
for f in brainfits:
 edges.append({'source':'assay:NIAGADS:ROSMAP_CUIMC1_2_MIT:SMAD5:exc','target':'study:Mak_FI','relation':'disease_molecular_colocalization_sensitivity','method':'SuSiE-coloc','nsnps':int(f['nsnps']),'p1':.0001,'p2':.0001,'p12':float(f['p12']),'PP_H3':float(f['PP.H3.abf']),'PP_H4':float(f['PP.H4.abf']),'ld_reference':'1000 Genomes European 503 participants, rank-deficient','input_file':brainfit.name,'input_sha256':hashlib.sha256(brainfit.read_bytes()).hexdigest(),'caveat':'Postselected cis-window model, LD mismatch and incomplete support; not proof of distinct or shared causal signal.'})
nodefile=p/'stage6_evidence_graph_nodes.jsonl';edgefile=p/'stage6_evidence_graph_edges.jsonl'
for file,vals in [(nodefile,nodes.values()),(edgefile,edges)]:
 with file.open('w') as f:
  for x in vals:f.write(json.dumps(x,sort_keys=True)+'\n')
assert all(x['source'] in nodes and x['target'] in nodes for x in edges)
meta={'input_sha256':sha,'candidate_count':len(variant_ids),'candidate_lists':{'weighted_15':len(rows),'shared_22':len(rows22),'overlap':len(ids22.intersection(r['variant'] for r in rows))},'additional_input_sha256':sha22,'nodes':len(nodes),'edges':len(edges),'edge_counts':{k:sum(e['relation']==k for e in edges) for k in sorted(set(e['relation'] for e in edges))},'warning':'Descriptive selected-source graph, not a trained ML model, an independent validation, a clinical risk score, or evidence of causality. No AlphaGenome predictions or derived features included.'}
(p/'stage6_evidence_graph_receipt.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))

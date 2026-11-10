import csv,collections,json,pathlib,urllib.parse
q=list(csv.DictReader(open('outputs/stage5_score_queue.csv')));pro=list(csv.DictReader(open('outputs/stage5_molecular_credible_set_provenance.csv')));meas=list(csv.DictReader(open('outputs/stage5_measured_molecular_evidence.csv')));avi={x['variant_GRCh38_ref_alt']:x for x in csv.DictReader(open('outputs/stage5_23_markers_atlas_avi.csv'))}
for x in csv.DictReader(open('outputs/stage5_avi_matched_control_comparison.csv')):avi[x['control_variant']]={'avi_phred_displayed_rounded':x['control_avi_phred_displayed'],'top_genomewide_snv_percent_displayed':x['control_top_genomewide_percent_displayed'],'source_url':x['source_url']}
for x in csv.DictReader(open('outputs/stage5_priority_credible_snv_atlas_avi.csv')):avi[x['variant_GRCh38_ref_alt']]=x
measure=collections.defaultdict(list)
for x in meas:measure[x['variant_key']].append(x)
prior=collections.defaultdict(list)
for x in pro:prior[x['variant']].append(x)
markers={x['variant']:x for x in csv.DictReader(open('outputs/shared22_eight_source_evidence_matrix.csv'))};extra=next(x for x in csv.DictReader(open('outputs/integrated_15_candidate_review.csv')) if x['variant']=='rs323509');observations=collections.defaultdict(list)
for p in ['outputs/stage5_initial_atlas_tissue_observations.csv','outputs/stage5_rs323509_atlas_initial_observations.csv']:
 for x in csv.DictReader(open(p)):
  if x.get('output_type')=='AVI':continue
  observations[x['variant_id']].append(x)
rows=[]
for x in q:
 key=x['variant_key'];a=avi.get(key,{});labels=x['input_labels'].split('|');refs=[r for label in labels for r in prior.get(label,[])];m=measure.get(key,[]);mark=next((l for l in labels if l.startswith('rs') and l in markers),None);trait=markers.get(mark,{});obs=observations.get(key,[])
 if mark:traittext='; '.join(f'{t}: p={trait[t+"_p"]}' for t in ['ukb_insomnia','atkins_FI','ye_FFS'])
 elif 'rs323509' in labels:traittext=f"UKB insomnia p={extra['UKB_insomnia_p']}; Atkins FI p={extra['FI_p']}; Ye FFS p={extra['FFS_p']}"
 else:traittext='not assigned; molecular credible-set membership does not establish GWAS association'
 rows.append({'variant_GRCh38_ref_alt':key,'input_classes':x['input_classes'],'input_labels':x['input_labels'],'candidate_gene_context':x['gene_context'],'molecular_credible_set_ids':'|'.join(sorted({v['credible_set_id'] for v in refs if v['credible_set_id']})),'molecular_pips':'|'.join(sorted({v['molecular_pip'] for v in refs if v['molecular_pip']})),'measured_QTL_record_count':len(m),'measured_QTL_gene_ids':'|'.join(sorted({v['gene_id'] for v in m if v['gene_id']})),'measured_QTL_tissues':'|'.join(sorted({v['tissue_or_cell'] for v in m if v['tissue_or_cell']})),'disease_association_context':traittext,'atlas_avi_phred_displayed':a.get('avi_phred_displayed_rounded',''),'atlas_top_genomewide_percent_displayed':a.get('top_genomewide_snv_percent_displayed',''),'atlas_web_url':a.get('source_url',''),'tissue_prediction_observation_count':len(obs),'gene_expression_splice_full_precision':'not obtained','chromatin_full_precision':'not obtained','validation_status':'pending full-precision output and tissue/allele/event match','status':'precomputed AVI read from website' if a else 'no AlphaGenome prediction obtained','limitations':'postselected; molecular QTL not proof of disease causality; no modeled individual risk; study samples can overlap'})
with open('outputs/stage5_preliminary_functional_atlas.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print('rows',len(rows),'AVI obtained',sum(bool(x['atlas_avi_phred_displayed']) for x in rows),'tissue observations',sum(bool(x['tissue_prediction_observation_count']) for x in rows),'unscored',sum(not bool(x['atlas_avi_phred_displayed']) for x in rows))

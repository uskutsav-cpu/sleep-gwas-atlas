import csv,collections,json,hashlib,pathlib
r=list(csv.DictReader(open('outputs/stage5_variant_input_manifest.csv')));norm={x['source_variant_id']:x for x in csv.DictReader(open('outputs/stage5_indel_normalization.csv'))};groups=collections.defaultdict(list)
for x in r:
 if x['status'] not in ('reference_verified_study_pair','reference_matches_GRCh38'):continue
 n=norm.get(x['input_label']);pos=n['normalized_position_GRCh38'] if n else x['position_GRCh38'];ref=n['normalized_ref'] if n else x['reference_allele'];alt=n['normalized_alt'] if n else x['alternate_allele']
 if n and n['normalization_status']!='normalized_within_window':continue
 groups[(x['chr'],int(pos),ref,alt)].append(x)
controls=pathlib.Path('outputs/stage5_matched_controls.csv')
if controls.exists():
 for x in csv.DictReader(controls.open()):
  if x['status_GRCh38']!='unique_liftover_rsid_pair_reference_verified':continue
  key=(x['chr_GRCh38'],int(x['pos_GRCh38']),x['ref_GRCh38'],x['alt_GRCh38'])
  groups[key].append({'input_label':'control_for_'+x['candidate'],'input_class':'matched_control','aliases':x['control_rsid'],'gene_context':'','status':x['status_GRCh38']})
rows=[]
for (chrom,pos,ref,alt),items in groups.items():
 labels=sorted({x['input_label'] for x in items});classes=sorted({x['input_class'] for x in items});genes=sorted({g for x in items for g in x['gene_context'].split('|') if g});rows.append({'variant_key':f'{chrom}:{pos}:{ref}>{alt}','chr':chrom,'position_GRCh38':pos,'reference_allele':ref,'alternate_allele':alt,'input_classes':'|'.join(classes),'input_labels':'|'.join(labels),'aliases':'|'.join(sorted({a for x in items for a in x['aliases'].split('|') if a})),'gene_context':'|'.join(genes),'status':'reference_checked_not_scored','note':'molecular credible-set membership does not imply GWAS causality; candidate coordinate and REF validated separately'})
priority={'rs323509':0,'rs12536335':1,'rs12671330':1,'rs7734045':2}
def rank(x):
 labs=x['input_labels'].split('|');p=min([priority.get(l,10) for l in labs]);return (0 if 'exploratory_marker' in x['input_classes'] else (1 if 'matched_control' in x['input_classes'] else 2),p,x['chr'],int(x['position_GRCh38']),x['alternate_allele'])
rows.sort(key=rank)
path=pathlib.Path('outputs/stage5_score_queue.csv')
with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print('unique reference-checked allele records',len(rows),'markers',sum('exploratory_marker' in x['input_classes'] for x in rows),'indels',sum(len(x['reference_allele'])!=len(x['alternate_allele']) for x in rows))

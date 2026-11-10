import csv,collections,pathlib
rows=list(csv.DictReader(open('outputs/shared22_eight_source_evidence_matrix.csv'))); ex=next(x for x in csv.DictReader(open('outputs/integrated_15_candidate_review.csv')) if x['variant']=='rs323509');rows.append({'variant':'rs323509','chr_GRCh37_assumed':ex['GRCh37_chr'],'position_GRCh37_assumed':ex['GRCh37_pos'],'nearest_gene_body_positional':ex['nearest_gene_body_positional'],'atkins_FI_p':ex['FI_p'],'ye_FFS_p':ex['FFS_p'],'ukb_insomnia_p':ex['UKB_insomnia_p'],'finngen_insomnia_p':ex['FinnGen_insomnia_p'],'mvp_insomnia_p':ex['MVP_insomnia_p'],'hfrs_p':ex['FinnGen_HFRS_p']})
qt=list(csv.DictReader(open('outputs/stage3_gtex_v10_candidate_qtls.csv')));bychr=collections.defaultdict(list)
for x in rows:bychr[x['chr_GRCh37_assumed']].append(x)
components=[]
for chrom,arr in bychr.items():
 arr.sort(key=lambda x:int(x['position_GRCh37_assumed'])); group=[]
 for x in arr:
  if group and int(x['position_GRCh37_assumed'])-int(group[-1]['position_GRCh37_assumed'])>1000000:components.append(group);group=[]
  group.append(x)
 if group:components.append(group)
result=[]
for g in components:
 for x in g:x['joint_triage_max_p']=max(float(x[k]) for k in ('atkins_FI_p','ye_FFS_p','finngen_insomnia_p'))
 lead=min(g,key=lambda x:x['joint_triage_max_p']);vs={x['variant'] for x in g};q=[x for x in qt if x['candidate'] in vs];qgenes=sorted(set(x['gene_symbol'] for x in q if x['gene_symbol']));
 row={'locus':'chr'+g[0]['chr_GRCh37_assumed']+':'+str(min(int(x['position_GRCh37_assumed']) for x in g))+'-'+str(max(int(x['position_GRCh37_assumed']) for x in g))+' (GRCh37 assumed)','markers':';'.join(x['variant'] for x in g),'triage_marker':lead['variant'],'atkins_FI_p':lead['atkins_FI_p'],'ye_FFS_p':lead['ye_FFS_p'],'finngen_insomnia_p':lead['finngen_insomnia_p'],'ukb_insomnia_p':lead['ukb_insomnia_p'],'mvp_insomnia_p':lead['mvp_insomnia_p'],'finngen_HFRS_p':lead['hfrs_p'],'mak_FI_p':lead.get('mak_FI_p','not available in 22-marker matrix'),'mak_FP_p':lead.get('mak_FP_p','not available in 22-marker matrix'),'max_fi_ffs_finn_insomnia_p':lead['joint_triage_max_p'],'positional_gene_labels':';'.join(sorted(set(x['nearest_gene_body_positional'] for x in g))),'gtex_v10_exact_marker_qtl_gene_labels':';'.join(qgenes),'gtex_v10_eqtl_rows':sum(x['qtl_type']=='Eqtl' for x in q),'gtex_v10_sqtl_rows':sum(x['qtl_type']=='Sqtl' for x in q),'PLACO_plus':'Unverified source and QC; see PLACO audit; not used to rank','status':'Exploratory postselected candidate-centered triage; not independent loci or confirmed gene'}
 result.append(row)
result.sort(key=lambda x:x['max_fi_ffs_finn_insomnia_p']);
for i,x in enumerate(result,1):x['triage_rank']=i
assert len(rows)==23 and len(result)==19
p=pathlib.Path('outputs/stage3_provisional_locus_triage.csv')
with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=['triage_rank']+[k for k in result[0] if k!='triage_rank']);w.writeheader();w.writerows(result)
print('23 candidates; 19 connected groups; top five',[(x['triage_marker'],round(x['max_fi_ffs_finn_insomnia_p'],6)) for x in result[:5]])

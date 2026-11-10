import csv,collections
r=list(csv.DictReader(open('outputs/stage3_abc_ld_proxy_regulatory_links.csv')));genes=collections.defaultdict(dict)
for f in ['outputs/stage4_mak_fi_magma_gene_results.csv','outputs/stage4_ye_ffs_magma_gene_results.csv']:
 for x in csv.DictReader(open(f)):
  if x['gene_symbol_v92']:genes[x['gene_symbol_v92']][x['GENE']]=x
groups=collections.defaultdict(list)
for x in r:groups[(x['candidate'],x['TargetGene'])].append(x)
rows=[]
for (c,s),group in sorted(groups.items()):
 matches=list(genes[s].values());valid=[g for g in matches if all(g['CHR']==x['proxy_chr_GRCh37'] and int(g['START'])<=int(x['TargetGeneTSS'])<=int(g['STOP']) for x in group)]
 status='unique match among tested Ensembl v92 coding genes with matching chromosome and TSS' if len(matches)==len(valid)==1 else 'unresolved: absent, ambiguous, or coordinate discordant in tested coding-gene subset'
 rows.append({'candidate':c,'abc_target_symbol':s,'candidate_gene_id_v92':valid[0]['GENE'] if len(matches)==len(valid)==1 else '','mapping_status':status,'n_matching_tested_ids':len(matches),'source_tss_grch37':'|'.join(sorted({x['proxy_chr_GRCh37']+':'+x['TargetGeneTSS'] for x in group})),'prediction_records':len(group),'interpretation':'identifier mapping only; ABC and LD do not establish disease causality'})
with open('outputs/stage3_abc_ld_proxy_gene_id_mapping.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print('candidate-target labels',len(rows),'mapped',sum(bool(x['candidate_gene_id_v92']) for x in rows))

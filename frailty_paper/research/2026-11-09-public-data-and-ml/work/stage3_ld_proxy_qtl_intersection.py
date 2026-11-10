import csv,collections,pathlib,json
ps=list(csv.DictReader(open('outputs/stage3_candidate_reference_ld_proxies.csv')));by=collections.defaultdict(list)
for x in ps:
 if x['proxy_pos_GRCh37']!=x['candidate_pos_GRCh37']:by[(x['candidate'],x['proxy_reference_id'])].append(x)
work=list(csv.DictReader(open('outputs/stage3_remaining_molecular_worklist.csv')));jobs=[]
for x in work:
 if x['candidate'] and pathlib.Path('work/stage3_remaining_qtl/'+x['candidate']+'.tsv').exists():jobs.append((x['candidate'],x['gene_id'],x['gene_symbol'],x['tissue_v10'],x['dataset_r7'],pathlib.Path('work/stage3_remaining_qtl/'+x['candidate']+'.tsv')))
jobs.append(('rs323509','ENSG00000251574','ENSG00000251574','Testis','QTD000340',pathlib.Path('work/stage3_QTD000340_rs323509_testis_event.tsv')))
rows=[];checks=[]
for candidate,gid,symbol,tissue,dataset,path in jobs:
 q=collections.defaultdict(list)
 for x in csv.DictReader(path.open(),delimiter='\t'):
  if x['gene_id']==gid and x['rsid']:q[x['rsid']].append(x)
 matched=ambiguous=0
 for (c,rs),refs in by.items():
  if c!=candidate or rs not in q:continue
  for ref in refs:
   al=set(ref['proxy_bim_A1_A2'].split('/'));v=[x for x in q[rs] if {x['ref'],x['alt']}==al and x['chromosome']==ref['candidate_chr_GRCh37']]
   if len(v)!=1:ambiguous+=1;continue
   x=v[0];matched+=1
   rows.append({'candidate':candidate,'proxy_reference_id':rs,'proxy_GRCh37_chr':ref['candidate_chr_GRCh37'],'proxy_GRCh37_pos':ref['proxy_pos_GRCh37'],'proxy_allele_pair':ref['proxy_bim_A1_A2'],'candidate_proxy_r2':ref['r2_genotype_dosage'],'gene_id':gid,'source_gene_symbol':symbol,'tissue':tissue,'qtl_catalogue_release':'r7','dataset':dataset,'molecular_trait_id':x['molecular_trait_id'],'qtl_variant_grch38':x['variant'],'qtl_p_unadjusted':x['pvalue'],'qtl_beta_alt':x['beta'],'qtl_se':x['se'],'note':'nominal selected-region QTL, joined by rsID and allele pair across assumed builds; no multiple-testing correction, shared causality or allele-direction alignment to marker'})
 checks.append({'candidate':candidate,'gene_id':gid,'dataset':dataset,'source_qtl_rows':sum(len(a) for a in q.values()),'matched_proxy_records':matched,'ambiguous_or_allele_mismatch_records':ambiguous})
with open('outputs/stage3_selected_molecular_qtl_ld_proxies.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
json.dump({'analysis':'intersection of already selected regional QTL summaries with high LD nonself proxies; not a comprehensive QTL screen','joins':'same rsID, chromosome, reference/alternate unordered allele pair; assumed GRCh37→GRCh38 correspondence, no independent per-variant liftover; exactly one QTL match required','checks':checks,'output_rows':len(rows)},open('outputs/stage3_selected_molecular_qtl_ld_proxy_audit.json','w'),indent=2)
print('output',len(rows),'per candidate',dict(collections.Counter(x['candidate'] for x in rows)));print('p<.05 unadjusted',sum(float(x['qtl_p_unadjusted'])<.05 for x in rows))

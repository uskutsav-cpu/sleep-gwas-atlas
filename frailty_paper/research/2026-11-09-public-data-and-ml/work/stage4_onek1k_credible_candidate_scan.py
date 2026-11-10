import csv,gzip,hashlib,json,pathlib,collections
m=list(csv.DictReader(open('outputs/shared22_eight_source_evidence_matrix.csv')));ex=next(x for x in csv.DictReader(open('outputs/integrated_15_candidate_review.csv')) if x['variant']=='rs323509');targets=[{'candidate':x['variant'],'chrom':x['chr_GRCh38'],'pos':x['position_GRCh38'],'a':x['score_counted_allele'],'b':x['other_allele']} for x in m];targets.append({'candidate':'rs323509','chrom':ex['GRCh37_chr'],'pos':ex['GRCh38_pos'],'a':ex['counted_allele'],'b':ex['other_allele']}); assert len(targets)==23
lookup={(x['chrom'],x['pos'],frozenset([x['a'],x['b']])):x for x in targets};assert len(lookup)==23
sources={608:'B_naive',609:'CD14_Mono',612:'CD4_Naive'};out=[];receipts=[]
for i,cell in sources.items():
 p=pathlib.Path(f'work/QTD000{i}.credible_sets.tsv.gz');n=0
 with gzip.open(p,'rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   n+=1;parts=row['variant'].split('_')
   if len(parts)!=4:continue
   c,pos,ref,alt=parts;tg=lookup.get((c.removeprefix('chr'),pos,frozenset((ref,alt))))
   if tg:out.append({'candidate':tg['candidate'],'candidate_alleles':tg['a']+'/'+tg['b'],'qtl_variant_id':row['variant'],'qtl_rsid':row['rsid'],'cell_population':cell,'dataset':f'QTD000{i}','gene_id':row['gene_id'],'molecular_trait_id':row['molecular_trait_id'],'credible_set_id':row['cs_id'],'credible_set_size':row['cs_size'],'pip':row['pip'],'p_value':row['pvalue'],'beta':row['beta'],'se':row['se'],'cs_min_r2':row['cs_min_r2'],'region':row['region'],'caveat':'Molecular eQTL credible-set membership only; no GWAS colocalization, no trait causal assignment; original variant build and alleles require verification'})
 receipts.append({'dataset':f'QTD000{i}','cell_population':cell,'source_url':f'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000038/QTD000{i}/QTD000{i}.credible_sets.tsv.gz','source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'total_rows':n})
cols=['candidate','candidate_alleles','qtl_variant_id','qtl_rsid','cell_population','dataset','gene_id','molecular_trait_id','credible_set_id','credible_set_size','pip','p_value','beta','se','cs_min_r2','region','caveat'];p=pathlib.Path('outputs/stage4_onek1k_candidate_credible_set_matches.csv')
with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(out)
unique={(x['candidate'],x['qtl_variant_id'],x['dataset'],x['molecular_trait_id'],x['credible_set_id']) for x in out};receipt={'release':'eQTL Catalogue r7','study':'QTS000038 OneK1K','inputs':receipts,'candidates_queried':23,'raw_matches':len(out),'unique_variant_phenotype_set_matches':len(unique),'matched_candidates':sorted({x['candidate'] for x in out}),'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'limits':'Only selected molecular credible sets, not all significant or tested cell QTL; candidates selected previously; exact GWAS build and genotype contrasts unverified; multiple rsIDs can duplicate same molecular variant and set; not an enrichment or GWAS colocalization.'};pathlib.Path('outputs/stage4_onek1k_candidate_credible_set_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

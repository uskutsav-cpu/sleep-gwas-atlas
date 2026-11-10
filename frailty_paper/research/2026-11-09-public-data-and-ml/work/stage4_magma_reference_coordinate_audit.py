import csv,collections,pathlib,json
sources={'Mak_FI':('/Volumes/Extreme SSD/Codex Archive/mak2023_same_cohort/fi.fastGWA','SNP','CHR','POS'),'Ye_FFS':('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv','SNP','CHROM','POS')}
want=collections.defaultdict(list);counts={}
for name,(path,snp,chrom,pos) in sources.items():
 with open(path) as f:
  for x in csv.DictReader(f,delimiter='\t'):
   rs=x[snp]
   if rs.startswith('rs') and rs[2:].isdigit() and int(rs[2:])%1000==0:want[rs].append((name,x[chrom],x[pos]))
ref=collections.defaultdict(list)
with open('work/stage4_magma_reference/g1000_eur.bim') as f:
 for l in f:
  x=l.split();rs=x[1]
  if rs in want:ref[rs].append((x[0],x[3]))
examples=collections.defaultdict(list)
for rs,uses in want.items():
 for name,ch,pos in uses:
  c=counts.setdefault(name,collections.Counter());c['sampled_source_rows']+=1
  if rs not in ref:c['not_directly_in_reference']+=1;continue
  if len(ref[rs])!=1:c['duplicate_reference_id']+=1;continue
  if (ch,pos)==ref[rs][0]:c['direct_id_and_coordinate_match']+=1
  else:
   c['direct_id_coordinate_mismatch']+=1
   if len(examples[name])<10:examples[name].append({'rsid':rs,'source':[ch,pos],'reference':ref[rs][0]})
out={'sampling_rule':'numeric rsID divisible by 1000; fixed non-random ID sample','results':{k:dict(v) for k,v in counts.items()},'mismatch_examples':examples,'scope':'direct rsID only; does not test synonyms, allele mapping, unsampled IDs or all ambiguity; source builds assumed GRCh37'}
pathlib.Path('outputs/stage4_magma_reference_coordinate_sample_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

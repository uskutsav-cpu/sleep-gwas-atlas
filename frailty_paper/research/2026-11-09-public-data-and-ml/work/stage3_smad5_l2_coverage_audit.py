import csv,gzip,pathlib,collections,json,hashlib,math
rows=list(csv.DictReader(open('outputs/stage3_skeletal_qtl_credible_set_variants.csv')));d={}
for x in rows:
 if x['qtl_phenotype']=='SMAD5 selected splice event' and x['cs_id'].endswith('_L2'):
  v=x['variant'];d.setdefault(v,{'variant':v,'rsids':set(),'pip':x['pip']})['rsids'].add(x['rsid'])
assert len(d)==18
smim=list(csv.DictReader(open('work/stage3_QTD000281_SMIM32_1mb.tsv'),delimiter='\t'));splice=list(csv.DictReader(open('work/stage3_QTD000285_SMAD5_splice_1mb.tsv'),delimiter='\t'));unique=collections.Counter(r['rsid'] for r in smim);smim_byvar=collections.defaultdict(list);splice_byvar=collections.defaultdict(list)
for r in smim:smim_byvar[r['variant']].append(r)
for r in splice:splice_byvar[r['variant']].append(r)
positions={int(v.split('_')[1]) for v in d};src=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw/finngen_R13_F5_INSOMNIA.gz');found=collections.defaultdict(list);n=0
with gzip.open(src,'rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  n+=1
  if x['#chrom']=='5' and int(x['pos']) in positions:found[int(x['pos'])].append(x)
print('scanned',n,'target positions returned',len(found),flush=True)
previous={r['variant'] for r in csv.DictReader(open('work/stage3_smad5_FinnGen_insomnia_aligned.csv'))};out=[]
for variant,a in d.items():
 chrom,pos,ref,alt=variant.split('_');snps=smim_byvar.get(variant,[]);sp=splice_byvar.get(variant,[]);target=[x for x in found[int(pos)] if x['ref']==ref and x['alt']==alt];pair={ref,alt};flags=[]
 if not (len(ref)==len(alt)==1 and pair<=set('ACGT')):flags.append('non_snp')
 elif pair in ({'A','T'},{'C','G'}):flags.append('palindromic')
 if not snps:flags.append('absent_sMIM32_expression_query')
 if snps and all(float(x['maf'])<.01 for x in snps):flags.append('expression_qtl_maf_lt_0.01')
 if snps and all(unique[x['rsid']]!=1 for x in snps):flags.append('expression_qtl_rsid_not_unique')
 if not target:flags.append('no_exact_finngen_position_ref_alt')
 elif len(target)>1:flags.append('duplicate_exact_finngen_rows')
 else:
  try:
   af=float(target[0]['af_alt'])
   if not .01<=af<=.99:flags.append('finngen_maf_lt_0.01')
  except:flags.append('finngen_frequency_missing_or_invalid')
 out.append({'variant_GRCh38':variant,'rsids':';'.join(sorted(a['rsids'])),'pip':a['pip'],'SMIM32_qtl_rows':len(snps),'SMIM32_qtl_mafs':';'.join(sorted({x['maf'] for x in snps})),'SMAD5_sqtl_rows':len(sp),'SMAD5_sqtl_mafs':';'.join(sorted({x['maf'] for x in sp})),'finngen_exact_rows':len(target),'finngen_af_alt':target[0]['af_alt'] if len(target)==1 else '', 'finngen_p':target[0]['pval'] if len(target)==1 else '', 'previous_input_matched_rsids':';'.join(sorted(a['rsids']&previous)),'observed_exclusion_flags':';'.join(flags)})
out.sort(key=lambda x:int(x['variant_GRCh38'].split('_')[1]));p=pathlib.Path('outputs/stage3_smad5_L2_finngen_coverage_audit.csv')
with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
rec={'variants':len(out),'source_scanned':n,'source_file':str(src),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'none_in_previous_input':all(not x['previous_input_matched_rsids'] for x in out),'flag_counts':dict(collections.Counter(y for x in out for y in x['observed_exclusion_flags'].split(';') if y)),'limitation':'Flags can co-occur; ID and variant representation differences may affect unique counts; no new colocalization or proof of lack of effect.'};pathlib.Path('outputs/stage3_smad5_L2_finngen_coverage_receipt.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2),flush=True)

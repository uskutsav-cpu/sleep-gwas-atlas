import gzip,csv,collections,hashlib,math,pathlib,json
p=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/v9-public-data/atkins-2021/Atkins_2021_FrailtyIndex_GWAS_EuropeansAged60to70years_GRCh37_UKBiobank_TwinGene_associations.txt.gz')
mat={x['variant']:x for x in csv.DictReader(open('outputs/six_candidates_six_trait_direction_matrix.csv'))};genes={x['variant']:x for x in csv.DictReader(open('outputs/song124_six_hit_positional_genes.csv'))};assert len(mat)==6
s=hashlib.sha256()
with p.open('rb') as f:
 for v in iter(lambda:f.read(8*1024*1024),b''):s.update(v)
found=collections.defaultdict(list);n=0
with gzip.open(p,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  n+=1
  if r['SNP'] in mat:found[r['SNP']].append(r)
result=[]
for rs,x in mat.items():
 g=genes[rs];hits=[y for y in found[rs] if y['CHR']==g['chr_GRCh37_assumed'] and y['BP']==g['position_GRCh37_assumed'] and {y['EA'],y['OA']}=={x['counted_allele'],x['other_allele']}]
 d={'variant':rs,'counted_allele':x['counted_allele'],'other_allele':x['other_allele'],'positional_gene':g['nearest_names'],'source_id_rows':len(found[rs]),'source_exact_rows':len(hits)}
 if len(hits)==1:
  h=hits[0];sign=1 if h['EA']==x['counted_allele'] else -1
  try:
   ub,us,tb,ts=sign*float(h['BETA_UKB1']),float(h['SE_UKB1']),-sign*float(h['BETA_TG']),float(h['SE_TG']);assert all(map(math.isfinite,(ub,us,tb,ts))) and us>0 and ts>0
   delta=tb-ub;se=math.hypot(us,ts);meta=(ub/us**2+tb/ts**2)/(1/us**2+1/ts**2);stored=float(x['FI_aligned_beta']);assert abs(meta-stored)<0.0001
   d.update(source_EA=h['EA'],source_OA=h['OA'],UKB_beta=ub,UKB_se=us,UKB_p=h['P_UKB1'],TwinGene_raw_beta=h['BETA_TG'],TwinGene_beta_adapter=tb,TwinGene_se=ts,TwinGene_p=h['P_TG'],TwinGene_minus_UKB=delta,difference_se_assuming_independence=se,difference_ci_low=delta-1.95996398454*se,difference_ci_high=delta+1.95996398454*se,difference_p_unadjusted=math.erfc(abs(delta/se)/math.sqrt(2)),reconstructed_meta_beta=meta,stored_meta_beta=stored)
  except (ValueError,AssertionError) as e:d['status']='invalid_values_or_reconstruction_mismatch'
 result.append(d)
out=pathlib.Path('outputs/six_leads_fi_cohort_comparison.csv');keys=list(dict.fromkeys(k for d in result for k in d))
with out.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(result)
receipt={'source_sha256':s.hexdigest(),'source_rows':n,'source_exact_unique_count':sum(x['source_exact_rows']==1 for x in result),'reconstructions_passed':sum('reconstructed_meta_beta' in x for x in result),'output_sha256':hashlib.sha256(out.read_bytes()).hexdigest()};pathlib.Path('outputs/six_leads_fi_cohort_comparison.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
for x in result:print(x)

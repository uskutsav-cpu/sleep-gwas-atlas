import gzip,csv,collections,hashlib,math,pathlib,json,pysam,numpy as np
src=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/v9-public-data/atkins-2021/Atkins_2021_FrailtyIndex_GWAS_EuropeansAged60to70years_GRCh37_UKBiobank_TwinGene_associations.txt.gz')
base='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
rows=list(csv.DictReader(open('outputs/expanded_22_candidate_risk_gene_table.csv'))); wanted={x['variant']:x for x in rows};found=collections.defaultdict(list);n=0
with gzip.open(src,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  n+=1
  if r['SNP'] in wanted:found[r['SNP']].append(r)
print('Scanned',n,'source rows',flush=True)
result=[]
for x in rows:
 hits=[h for h in found[x['variant']] if h['CHR']==x['GRCh37_chr_assumed'] and h['BP']==x['GRCh37_pos_assumed'] and {h['EA'],h['OA']}=={x['counted_allele'],x['other_allele']}]
 assert len(hits)==1,(x['variant'],len(hits))
 h=hits[0];sign=1 if h['EA']==x['counted_allele'] else -1
 ub,us,tb,ts=sign*float(h['BETA_UKB1']),float(h['SE_UKB1']),-sign*float(h['BETA_TG']),float(h['SE_TG'])
 assert all(math.isfinite(v) for v in (ub,us,tb,ts)) and us>0 and ts>0
 meta=(ub/us**2+tb/ts**2)/(1/us**2+1/ts**2);stored=float(x['FI_beta_aligned']);assert abs(meta-stored)<0.0001,(x['variant'],meta,stored)
 result.append({'variant':x['variant'],'counted_allele':x['counted_allele'],'other_allele':x['other_allele'],'ukb_beta':ub,'ukb_se':us,'twingene_beta_inferred':tb,'twingene_se':ts,'meta_beta':stored,'reconstructed_meta_beta':meta,'ukb_vs_twingene_opposite_sign':ub*tb<0})
out=pathlib.Path('outputs/22_fi_cohort_weight_sensitivity.csv')
with out.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=result[0]);w.writeheader();w.writerows(result)
panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t'));eur={r['sample'] for r in panel if r['super_pop']=='EUR'};files={};data=[];calls=[]
for x in rows:
 chrom=x['GRCh37_chr_assumed'];pos=int(x['GRCh37_pos_assumed']);a=x['counted_allele'];b=x['other_allele']
 if chrom not in files:files[chrom]=pysam.VariantFile(base.format(chrom=chrom))
 f=files[chrom];hits=[v for v in f.fetch(chrom,pos-1,pos) if v.id==x['variant'] and v.pos==pos and len(v.alleles)==2 and set(v.alleles)=={a,b}];assert len(hits)==1
 v=hits[0];idx=v.alleles.index(a);d={}
 for s in sorted(eur & set(f.header.samples)):
  gt=v.samples[s].get('GT')
  if gt and len(gt)==2 and all(g is not None for g in gt):d[s]=sum(g==idx for g in gt)
 data.append(d);calls.append(len(d))
samples=sorted(set.intersection(*(set(d) for d in data)));assert len(samples)==503
M=np.array([[d[s] for d in data] for s in samples],float);W=np.array([[x['ukb_beta'],x['meta_beta'],x['twingene_beta_inferred']] for x in result]);S=M@W
res={'source_rows_scanned':n,'source_exact_matches':len(result),'n_reference_eur':len(samples),'score_order':['UKB_FI','meta_FI','TwinGene_FI_inferred'],'score_correlations':np.corrcoef(S,rowvar=False).tolist(),'score_sd':S.std(axis=0,ddof=1).tolist(),'n_opposite_ukb_twingene':sum(r['ukb_vs_twingene_opposite_sign'] for r in result),'variant_calls':calls,'source_csv_sha256':hashlib.sha256(pathlib.Path('outputs/expanded_22_candidate_risk_gene_table.csv').read_bytes()).hexdigest(),'sensitivity_csv_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'warning':'Selected marginal weights; genomewide sign adapter inferred, 503 potentially related reference individuals, no phenotypes; cohort correlations not predictive validation.'}
pathlib.Path('outputs/22_fi_cohort_weight_sensitivity.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))

import csv,pysam,numpy as np,json,hashlib
base='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
rows=list(csv.DictReader(open('outputs/unvalidated_15_variant_score_prototype.csv')))
with open('work/1kg_phase3_panel.tsv') as f:panel=list(csv.DictReader(f,delimiter='\t'))
eurs={r['sample'] for r in panel if r['super_pop']=='EUR'}; bypop={r['sample']:r['pop'] for r in panel}; files={};genotypes={};check=[]
for row in rows:
 rs=row['variant']; chrom=row['chr']; pos=int(row['position_GRCh37_inferred']);a=row['counted_allele'];b=row['other_allele']
 try:
  if chrom not in files: files[chrom]=pysam.VariantFile(base.format(chrom=chrom))
  f=files[chrom]; samples=sorted(eurs&set(f.header.samples));match=[x for x in f.fetch(chrom,pos-1,pos) if x.id==rs and x.pos==pos]
  if len(match)!=1:raise ValueError('matching records '+str(len(match)))
  x=match[0]; alleles=x.alleles
  if len(alleles)!=2 or set(alleles)!={a,b}: raise ValueError('biallelic contrast mismatch '+str(alleles))
  idx=alleles.index(a);calls={}
  for s in samples:
   gt=x.samples[s].get('GT')
   if gt and len(gt)==2 and all(g is not None for g in gt):calls[s]=sum(g==idx for g in gt)
  genotypes[rs]=calls;check.append({'variant':rs,'status':'matched','alleles':list(alleles),'n_called':len(calls)})
  print(rs,'matched',len(calls),flush=True)
 except Exception as e:check.append({'variant':rs,'status':str(e)});print(rs,'error',str(e),flush=True)
if len(genotypes)!=len(rows):
 json.dump({'status':'incomplete_no_scores_calculated','checks':check},open('outputs/reference_15_score_diagnostic.json','w'),indent=2);raise SystemExit(1)
samples=sorted(set.intersection(*(set(x) for x in genotypes.values())))
M=np.array([[genotypes[r['variant']][s] for r in rows] for s in samples],dtype=float)
cols=['FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight'];weights=np.array([[float(r[k]) for k in cols] for r in rows]);scores=M@weights
loo=[]
for j,row in enumerate(rows):
 z=scores-np.outer(M[:,j],weights[j,:]); q=np.corrcoef(z,rowvar=False)
 loo.append({'omitted_variant':row['variant'],'FI_FFS_r':float(q[0,1]),'FI_insomnia_r':float(q[0,2]),'FFS_insomnia_r':float(q[1,2])})
subgroups={}
for pop in sorted({bypop[s] for s in samples}):
 ids=[i for i,s in enumerate(samples) if bypop[s]==pop]; mat=np.corrcoef(scores[ids,:],rowvar=False)
 subgroups[pop]={'n':len(ids),'correlations':mat.tolist()}
result={'subgroup_correlations':subgroups,'leave_one_out_correlations':loo,'status':'exploratory_reference_genotype_only_no_outcomes','sample_count':len(samples),'population_counts':{p:sum(bypop[s]==p for s in samples) for p in sorted({bypop[s] for s in samples})},'variant_count':len(rows),'checks':check,'score_columns':cols,'score_means':scores.mean(axis=0).tolist(),'score_sd_sample':scores.std(axis=0,ddof=1).tolist(),'score_pearson_correlations':np.corrcoef(scores,rowvar=False).tolist(),'source_score_sha256':hashlib.sha256(open('outputs/unvalidated_15_variant_score_prototype.csv','rb').read()).hexdigest(),'limitations':'This reference panel has no FI, FFS or insomnia outcomes. Scores are unvalidated and constructed from selected GWAS effects; no individual-level genotypes or scores are saved. Signed SNP effects are inferred harmonizations. EUR population structure can influence correlations.'}
print(json.dumps(result,indent=2),flush=True)
with open('outputs/reference_15_score_diagnostic.json','w') as f:json.dump(result,f,indent=2)

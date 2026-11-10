import csv,pysam,numpy as np,json,hashlib,pathlib
base='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
rows=list(csv.DictReader(open('outputs/unvalidated_15_variant_score_prototype.csv')))
panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t')); eur={r['sample'] for r in panel if r['super_pop']=='EUR'}; pop={r['sample']:r['pop'] for r in panel}
files={};geno={}
for row in rows:
 chrom=row['chr']; pos=int(row['position_GRCh37_inferred']); a=row['counted_allele'];b=row['other_allele'];rs=row['variant']
 if chrom not in files:files[chrom]=pysam.VariantFile(base.format(chrom=chrom))
 f=files[chrom]; found=[x for x in f.fetch(chrom,pos-1,pos) if x.id==rs and x.pos==pos];assert len(found)==1
 x=found[0];assert len(x.alleles)==2 and set(x.alleles)=={a,b};idx=x.alleles.index(a);calls={}
 for s in sorted(eur & set(f.header.samples)):
  gt=x.samples[s].get('GT')
  if gt and len(gt)==2 and all(g is not None for g in gt): calls[s]=sum(g==idx for g in gt)
 geno[rs]=calls; print(rs,len(calls),flush=True)
samples=sorted(set.intersection(*(set(x) for x in geno.values())));assert len(samples)==503
M=np.array([[geno[r['variant']][s] for r in rows] for s in samples],float)
cols=['FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight'];W=np.array([[float(r[k]) for k in cols] for r in rows]);S=M@W; old=json.load(open('outputs/reference_15_score_diagnostic.json')); assert np.allclose(np.corrcoef(S,rowvar=False),old['score_pearson_correlations'],rtol=0,atol=1e-12)

def summarize(X):
 C=np.cov(X,rowvar=False,ddof=1);result={}
 for j,name in enumerate(('FI','FFS','insomnia')):
  w=W[:,j]; V=float(w@C@w); diag=(w*w*np.diag(C))/V; contrib=(w*(C@w))/V
  assert np.isclose(sum(contrib),1);assert np.isclose(V,np.var(X@w,ddof=1));result[name]={'score_variance':V,'diagonal_fraction':diag.tolist(),'allocated_covariance_fraction':contrib.tolist(),'off_diagonal_fraction':float(1-sum(diag))}
 return result
pooled=summarize(M);sub={p:summarize(M[[i for i,s in enumerate(samples) if pop[s]==p],:]) for p in sorted(set(pop[s] for s in samples))}
rec=[]
for i,row in enumerate(rows):
 x={'variant':row['variant'],'counted_allele':row['counted_allele'],'reference_counted_allele_frequency':float(M[:,i].mean()/2)}
 for name in pooled:
  x[name+'_standalone_variance_fraction']=pooled[name]['diagonal_fraction'][i];x[name+'_allocated_covariance_fraction']=pooled[name]['allocated_covariance_fraction'][i]
 rec.append(x)
out=pathlib.Path('outputs/reference_15_score_variance_decomposition.csv')
with out.open('w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=rec[0]);wr.writeheader();wr.writerows(rec)
result={'n':len(samples),'input_score_sha256':hashlib.sha256(pathlib.Path('outputs/unvalidated_15_variant_score_prototype.csv').read_bytes()).hexdigest(),'output_csv_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'variants':[r['variant'] for r in rows],'pooled':pooled,'subgroups':sub,'note':'Descriptive variation in selected genotype sums; no outcomes; allocate covariance symmetrically, can be negative; samples possibly related, LD and population structure; not variance explained in disease'}
pathlib.Path('outputs/reference_15_score_variance_decomposition.json').write_text(json.dumps(result,indent=2)+'\n')
for name in pooled:
 a=sorted(rec,key=lambda x:x[name+'_allocated_covariance_fraction'],reverse=True);print(name,'offdiagonal',pooled[name]['off_diagonal_fraction'],'top',[(x['variant'],round(x[name+'_allocated_covariance_fraction'],4),round(x[name+'_standalone_variance_fraction'],4)) for x in a[:5]])

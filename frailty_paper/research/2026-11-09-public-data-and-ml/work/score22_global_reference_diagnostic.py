import csv,pathlib,pysam,numpy as np,json,hashlib
rows=list(csv.DictReader(open('outputs/unvalidated_22_variant_score_prototype.csv')));panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t'));pop={r['sample']:r['super_pop'] for r in panel};allids=set(pop);base='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';files={};geno=[];checks=[]
for r in rows:
 c=r['chr'];p=int(r['position_GRCh37_assumed']);a=r['counted_allele'];b=r['other_allele']
 if c not in files:files[c]=pysam.VariantFile(base.format(chrom=c))
 f=files[c];v=[x for x in f.fetch(c,p-1,p) if x.id==r['variant'] and x.pos==p and len(x.alleles)==2 and set(x.alleles)=={a,b}];assert len(v)==1
 x=v[0];ix=x.alleles.index(a);d={}
 for s in allids & set(f.header.samples):
  gt=x.samples[s].get('GT')
  if gt and len(gt)==2 and all(g is not None for g in gt):d[s]=sum(g==ix for g in gt)
 geno.append(d);checks.append({'variant':r['variant'],'calls':len(d)})
samples=sorted(set.intersection(*(set(g) for g in geno)));X=np.array([[g[s] for g in geno] for s in samples],float);names=('FI','FFS','insomnia');cols=('FI_beta_weight','FFS_beta_weight_inferred','UKB_insomnia_beta_weight');W=np.array([[float(r[c]) for c in cols] for r in rows]);scores=X@W;groups={}
for p in sorted(set(pop[s] for s in samples)):
 ix=[i for i,s in enumerate(samples) if pop[s]==p];group=X[ix];z=scores[ix];groups[p]={'n':len(ix),'score_sd':{name:float(np.std(z[:,j],ddof=1)) for j,name in enumerate(names)},'score_correlations':np.corrcoef(z,rowvar=False).tolist(),'counted_allele_frequencies':(group.mean(axis=0)/2).tolist()}
old=np.array(json.load(open('outputs/reference_22_vs_15_score_diagnostic.json'))['new22_cross_trait_score_correlations']);assert groups['EUR']['n']==503 and np.allclose(groups['EUR']['score_correlations'],old,rtol=0,atol=1e-12)
out=[]
for i,r in enumerate(rows):
 f={p:groups[p]['counted_allele_frequencies'][i] for p in groups};out.append({'variant':r['variant'],'counted_allele':r['counted_allele'],**{'reference_'+p+'_frequency':v for p,v in f.items()},'largest_minus_smallest_frequency':max(f.values())-min(f.values())})
p=pathlib.Path('outputs/22_score_reference_group_allele_frequencies.csv')
with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
result={'n_reference_samples_complete':len(samples),'groups':groups,'variant_order':[r['variant'] for r in rows],'per_variant_calls':checks,'weights_sha256':hashlib.sha256(pathlib.Path('outputs/unvalidated_22_variant_score_prototype.csv').read_bytes()).hexdigest(),'frequency_csv_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'notes':'1000 Genomes labels are broad heterogeneous convenience samples, not biological categories; 22 post-selected variants and European GWAS weights; no outcomes, predictive portability, calibration, clinical risks or causality inferred; small group sizes, admixture, potential relatives, and exact study-specific variant coding remain limitations. No individual genotypes/scores retained.'};pathlib.Path('outputs/22_score_reference_group_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n');print('n',len(samples));
for p,g in groups.items():print(p,g['n'],'SD',g['score_sd'],'correlation',g['score_correlations'][0][1])
for x in sorted(out,key=lambda x:x['largest_minus_smallest_frequency'],reverse=True)[:4]:print(x['variant'],round(x['largest_minus_smallest_frequency'],3))

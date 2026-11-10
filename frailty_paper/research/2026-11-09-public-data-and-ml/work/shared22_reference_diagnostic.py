import csv,pathlib,pysam,numpy as np,json,hashlib
shared=list(csv.DictReader(open('outputs/unvalidated_22_shared_direction_score.csv')));weighted={r['variant']:r for r in csv.DictReader(open('outputs/unvalidated_22_variant_score_prototype.csv'))};panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t'));eur={r['sample'] for r in panel if r['super_pop']=='EUR'};base='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';files={};data=[]
for r in shared:
 c=r['chr_GRCh37_assumed'];p=int(r['position_GRCh37_assumed']);a=r['counted_allele'];b=r['other_allele'];w=weighted[r['variant']];assert {a,b}=={w['counted_allele'],w['other_allele']}
 if c not in files:files[c]=pysam.VariantFile(base.format(chrom=c))
 f=files[c];hits=[v for v in f.fetch(c,p-1,p) if v.id==r['variant'] and v.pos==p and len(v.alleles)==2 and set(v.alleles)=={a,b}];assert len(hits)==1
 v=hits[0];ix=v.alleles.index(a);d={}
 for s in eur&set(f.header.samples):
  gt=v.samples[s].get('GT')
  if gt and len(gt)==2 and all(g is not None for g in gt):d[s]=sum(g==ix for g in gt)
 data.append(d)
samples=sorted(set.intersection(*(set(d) for d in data)));assert len(samples)==503
X=np.array([[d[s] for d in data] for s in samples],float);unweighted=X.sum(axis=1);W=[];G=[]
for i,r in enumerate(shared):
 w=weighted[r['variant']];oriented=X[:,i] if w['counted_allele']==r['counted_allele'] else 2-X[:,i];G.append(oriented);W.append([float(w[c]) for c in ('FI_beta_weight','FFS_beta_weight_inferred','UKB_insomnia_beta_weight')])
A=np.column_stack(G)@np.array(W);prior=json.load(open('outputs/reference_22_vs_15_score_diagnostic.json'));assert np.allclose(np.corrcoef(A,rowvar=False),prior['new22_cross_trait_score_correlations'],rtol=0,atol=1e-12)
corr=np.corrcoef(np.column_stack((unweighted,A)),rowvar=False);r={'n_reference_eur':503,'counted_variants':22,'shared_allele_count_mean':float(unweighted.mean()),'shared_allele_count_sd':float(unweighted.std(ddof=1)),'correlations_with_trait_sums':{'FI':float(corr[0,1]),'FFS':float(corr[0,2]),'insomnia':float(corr[0,3])},'shared_score_sha256':hashlib.sha256(pathlib.Path('outputs/unvalidated_22_shared_direction_score.csv').read_bytes()).hexdigest(),'caveat':'Genotype-only arithmetic in selected reference samples, no outcomes, LD/selection/encoding issues, no clinical risk or predictive validity. No individual genotypes or scores retained.'};pathlib.Path('outputs/shared22_reference_diagnostic.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

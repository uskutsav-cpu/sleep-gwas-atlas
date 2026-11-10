import csv,pathlib,pysam,numpy as np,json,hashlib,itertools
score=pathlib.Path('outputs/unvalidated_22_shared_direction_score.csv')
rows=list(csv.DictReader(score.open())); panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t'))
eur={x['sample'] for x in panel if x['super_pop']=='EUR'}
base='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
files={}; data=[]
for row in rows:
 c=row['chr_GRCh37_assumed']; p=int(row['position_GRCh37_assumed']); a=row['counted_allele']; b=row['other_allele']
 if c not in files: files[c]=pysam.VariantFile(base.format(chrom=c))
 f=files[c]; matches=[v for v in f.fetch(c,p-1,p) if v.id==row['variant'] and v.pos==p and len(v.alleles)==2 and set(v.alleles)=={a,b}]; assert len(matches)==1
 v=matches[0]; ix=v.alleles.index(a); d={}
 for s in eur & set(f.header.samples):
  gt=v.samples[s].get('GT')
  if gt and len(gt)==2 and all(k is not None for k in gt): d[s]=sum(k==ix for k in gt)
 data.append(d)
 print(row['variant'],flush=True)
samples=sorted(set.intersection(*(set(d) for d in data))); assert len(samples)==503
X=np.array([[d[s] for d in data] for s in samples],float); C=np.cov(X,rowvar=False,ddof=1); var=float(np.var(X.sum(axis=1),ddof=1)); assert np.isclose(C.sum(),var)
pairs=[]
for i,j in itertools.combinations(range(len(rows)),2):
 x,y=rows[i],rows[j]; same=x['chr_GRCh37_assumed']==y['chr_GRCh37_assumed']; distance=abs(int(x['position_GRCh37_assumed'])-int(y['position_GRCh37_assumed'])) if same else None
 pairs.append(dict(variant_1=x['variant'],variant_2=y['variant'],same_chromosome=same,distance_bp=distance,within_1mb=same and distance<=1000000,dosage_covariance=float(C[i,j]),variance_contribution=float(2*C[i,j])))
near=sum(p['variance_contribution'] for p in pairs if p['within_1mb']); far=sum(p['variance_contribution'] for p in pairs if not p['within_1mb']); diag=float(np.trace(C)); assert np.isclose(diag+near+far,var)
old=json.load(open('outputs/shared22_reference_diagnostic.json')); assert np.isclose(var,old['shared_allele_count_sd']**2,atol=1e-12)
r={'n_samples':len(samples),'n_markers':len(rows),'score_sha256':hashlib.sha256(score.read_bytes()).hexdigest(),'variance':var,'diagonal':diag,'within_1mb_covariance_contribution':near,'other_pairs_covariance_contribution':far,'within_1mb_fraction_of_total':near/var,'other_pairs_fraction_of_total':far/var,'pairs':pairs,'notes':'Sample dosage covariance, not phased haplotype r; descriptive, pooled EUR reference, no outcomes; no genotype records retained.'}
pathlib.Path('outputs/shared22_count_variance_diagnostic.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({k:v for k,v in r.items() if k!='pairs'},indent=2))

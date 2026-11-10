"""Genotype-only panel comparison; not outcome validation or chip imputation assessment."""
import csv, json, hashlib
from pathlib import Path
import numpy as np, pysam
O=Path('outputs')
old=list(csv.DictReader((O/'unvalidated_15_variant_score_prototype.csv').open()))
count=list(csv.DictReader((O/'unvalidated_22_shared_direction_score.csv').open()))
weighted=list(csv.DictReader((O/'unvalidated_22_variant_score_prototype.csv').open()))
chip=list(csv.DictReader((O/'brhs_metabochip_score_marker_screen.csv').open()))
panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t'))
samples=sorted(r['sample'] for r in panel if r['super_pop']=='EUR')
base='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
files={}; data={}; checks=[]
for rows,chrkey,poskey in [(old,'chr','position_GRCh37_inferred'),(count,'chr_GRCh37_assumed','position_GRCh37_assumed')]:
 for row in rows:
  v=row['variant']; a=row['counted_allele']; b=row['other_allele']; key=(v,a,b)
  if key in data: continue
  chrom=row[chrkey];pos=int(row[poskey]);
  if chrom not in files: files[chrom]=pysam.VariantFile(base.format(chrom=chrom))
  f=files[chrom]; hits=[x for x in f.fetch(chrom,pos-1,pos) if x.id==v and x.pos==pos and len(x.alleles)==2 and set(x.alleles)=={a,b}]
  assert len(hits)==1,(v,len(hits)); x=hits[0]; ix=x.alleles.index(a)
  arr=[]
  for sample in samples:
   gt=x.samples[sample].get('GT'); assert gt is not None and len(gt)==2 and all(g is not None for g in gt)
   arr.append(sum(g==ix for g in gt))
  data[key]=np.asarray(arr,float); checks.append(v)

def dosage(row):
 key=(row['variant'],row['counted_allele'],row['other_allele'])
 if key in data:return data[key]
 rev=(row['variant'],row['other_allele'],row['counted_allele'])
 assert rev in data
 return 2-data[rev]

def score(rows,field, subset=None):
 return sum((dosage(r)*float(r[field]) for r in rows if subset is None or r['variant'] in subset),np.zeros(len(samples)))
# Check that the regenerated reference sums reproduce the previous full-score comparisons.
prior=json.loads((O/'reference_22_vs_15_score_diagnostic.json').read_text())['same_trait_correlations_old15_vs_new22']
for name, a, b in [('FI','FI_beta_weight','FI_beta_weight'),('FFS','FFS_beta_weight_inferred','FFS_beta_weight_inferred'),('UKB_insomnia','insomnia_beta_weight','UKB_insomnia_beta_weight')]:
 observed=float(np.corrcoef(score(old,a),score(weighted,b))[0,1])
 assert np.isclose(observed,prior[name],rtol=0,atol=1e-12),(name,observed,prior[name])
results={}
for label,rows,fields in [('15',old,['FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight']),('22count',count,['weight']),('22weighted',weighted,['FI_beta_weight','FFS_beta_weight_inferred','UKB_insomnia_beta_weight'])]:
 cid='15' if label=='15' else '22'
 subset={r['variant'] for r in chip if r['score']==cid and r['build_matched_coordinate_match']=='True'}
 results[label]={'retained':sorted(subset),'weights':{}}
 for field in fields:
  full=score(rows,field); part=score(rows,field,subset)
  results[label]['weights'][field]={'pearson_r':float(np.corrcoef(full,part)[0,1]),'partial_over_full_sd':float(np.std(part,ddof=1)/np.std(full,ddof=1))}
res={'n_eur':len(samples),'reference':'1000 Genomes phase 3, EUR-labelled, pooled','results':results,'meaning':'Genotype-only reference comparison, not participant outcomes, imputation coverage, prediction, or evidence of marker quality.', 'chip_csv_sha256':hashlib.sha256((O/'brhs_metabochip_score_marker_screen.csv').read_bytes()).hexdigest(), 'score_sha256':{str(p):hashlib.sha256((O/p).read_bytes()).hexdigest() for p in ['unvalidated_15_variant_score_prototype.csv','unvalidated_22_shared_direction_score.csv','unvalidated_22_variant_score_prototype.csv']}}
(O/'brhs_chip_subset_reference.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))

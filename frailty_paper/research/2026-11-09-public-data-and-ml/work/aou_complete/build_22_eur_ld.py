import csv,json,hashlib
from pathlib import Path
import numpy as np,pysam
f=Path('outputs/unvalidated_22_shared_direction_score.csv')
r=list(csv.DictReader(f.open()));assert len(r)==22
panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t'));ids=sorted(x['sample'] for x in panel if x['super_pop']=='EUR');assert len(ids)==503
base='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';handles={};columns=[]
for x in r:
 c=x['chr_GRCh37_assumed'];p=int(x['position_GRCh37_assumed']);
 if c not in handles:handles[c]=pysam.VariantFile(base.format(c=c))
 hits=[v for v in handles[c].fetch(c,p-1,p) if v.id==x['variant'] and v.pos==p and len(v.alleles)==2 and set(v.alleles)=={x['counted_allele'],x['other_allele']}];assert len(hits)==1,(x['variant'],len(hits))
 rec=hits[0];a=rec.alleles.index(x['counted_allele']);vals=[]
 for person in ids:
  gt=rec.samples[person].get('GT');assert gt is not None and len(gt)==2 and all(y is not None for y in gt)
  vals.append(sum(y==a for y in gt))
 columns.append(vals);print(x['variant'],flush=True)
X=np.array(columns,float).T;R=np.corrcoef(X,rowvar=False);chrom=np.array([x['chr_GRCh37_assumed'] for x in r]);R=np.where(chrom[:,None]==chrom[None,:],R,0.)
assert np.allclose(R,R.T) and np.linalg.eigvalsh(R).min()>-1e-8
p=Path('outputs/stage6_aou_count22_EUR503_ld.csv')
with p.open('w',newline='') as h:
 w=csv.writer(h);w.writerow(['variant']+[x['variant'] for x in r]);w.writerows([[x['variant'],*row] for x,row in zip(r,R)])
Path('outputs/stage6_aou_count22_EUR503_ld_receipt.json').write_text(json.dumps({'n':503,'variants':22,'score_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'reference':'1000 Genomes phase 3 EUR, exact biallelic listed contrast and called genotypes; cross-chromosomal sample correlations set to zero','limitations':'May contain relatives and structure; remote full VCF checksum not independently verified; reference LD differs from AoU'},indent=2)+'\n')

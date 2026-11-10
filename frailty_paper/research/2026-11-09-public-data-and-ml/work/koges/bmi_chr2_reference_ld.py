import csv,json,hashlib
from pathlib import Path
import numpy as np,pysam
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
panel=list(csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t'))
groups={x['sample']:x['pop'] for x in panel if x['super_pop']=='EAS'}
vcf=pysam.VariantFile(url); samples=sorted(set(groups)&set(vcf.header.samples))
assert len(samples)==504
variants=[('rs35881094',58922921,'G','T'),('rs17550821',58978523,'A','G')]
mat=[]; metadata=[]
for rs,pos,counted,other in variants:
 matches=[v for v in vcf.fetch('2',pos-1,pos) if v.id==rs and v.pos==pos]
 if len(matches)!=1:raise RuntimeError('ambiguous variant '+rs)
 v=matches[0]
 if set(v.alleles)!={counted,other} or len(v.alleles)!=2:raise RuntimeError('variant allele mismatch '+rs)
 ix=v.alleles.index(counted); data=[]
 for s in samples:
  call=v.samples[s];gt=call.get('GT')
  if not gt or len(gt)!=2 or any(i is None for i in gt) or not call.phased:raise RuntimeError('missing/unphased '+rs+' '+s)
  data.append([int(i==ix) for i in gt])
 a=np.array(data);mat.append(a);metadata.append({'rsid':rs,'position_GRCh37':pos,'reference':v.ref,'alternate':','.join(v.alts),'counted':counted,'counted_allele_frequency':float(a.mean())})
def calc(indices):
 a=mat[0][indices].flatten();b=mat[1][indices].flatten()
 if np.std(a)==0 or np.std(b)==0:return {'n':len(indices),'status':'zero variance'}
 r=float(np.corrcoef(a,b)[0,1]);return {'n':len(indices),'r':r,'r2':r*r,'af_variant1':float(a.mean()),'af_variant2':float(b.mean()),'status':'ok'}
allidx=np.arange(len(samples)); main=calc(allidx)
rng=np.random.default_rng(20261029);boot=[];invalid=0
for i in range(2000):
 x=calc(rng.choice(allidx,len(allidx),replace=True))
 if x['status']=='ok':boot.append(x['r2'])
 else:invalid+=1
out={'source':url,'panel_sha256':hashlib.sha256(Path('work/1kg_phase3_panel.tsv').read_bytes()).hexdigest(),'reference':'1000 Genomes phase 3 East Asian superpopulation, not a Korean-specific panel','alleles':metadata,'pooled':main,'subpopulations':{pop:calc(np.array([i for i,s in enumerate(samples) if groups[s]==pop])) for pop in sorted(set(groups.values()))},'bootstrap_2000_individual_resamples_r2_quantiles_2_5_97_5':np.quantile(boot,[.025,.975]).tolist(),'bootstrap_invalid':invalid,'limitations':'Exploratory external reference, not KoGES LD or conditional association; subgroup estimates unstable, especially at low minor allele counts.'}
Path('outputs/stage6_koges_bmi_chr2_eas_reference_ld.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

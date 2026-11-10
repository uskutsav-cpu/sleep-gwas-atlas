#!/usr/bin/env python3
"""Exploratory matched-marker summary-Z sensitivity; not target score regression."""
import csv,hashlib,json,math
from pathlib import Path
from statistics import NormalDist
import numpy as np
import pysam
root=Path(__file__).resolve().parent.parent;out=root/'outputs'
def read(name,sep=','):
 with (out/name).open(newline='') as f:return list(csv.DictReader(f,delimiter=sep))
score=read('unvalidated_22_shared_direction_score.csv')
bbj={r['variant']:r for r in read('stage6_bbj_insomnia_oriented_candidates.csv') if r['score']=='shared22' and r['status']=='exact position and pair'}
prior=json.loads((out/'stage6_cross_population_matched14.json').read_text());selected=[r for r in score if r['variant'] in prior['markers']];assert [r['variant'] for r in selected]==prior['markers'] and len(selected)==14
panel=[]
with (root/'work/1kg_phase3_panel.tsv').open(newline='') as f:panel=list(csv.DictReader(f,delimiter='\t'))
names=sorted(r['sample'] for r in panel if r['super_pop']=='EAS');jpt={r['sample'] for r in panel if r['pop']=='JPT'};mask=np.array([n in jpt for n in names]);assert len(names)==504 and sum(mask)==104
handles={};vectors=[];components=[]
for r in selected:
 vname=r['variant'];b=bbj[vname];chrom=r['chr_GRCh37_assumed'];pos=int(r['position_GRCh37_assumed']);count=r['counted_allele'];other=r['other_allele'];assert (b['counted_allele'],b['other_allele'])==(count,other)
 if chrom not in handles:handles[chrom]=pysam.VariantFile(f'https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz')
 variants=[v for v in handles[chrom].fetch(chrom,pos-1,pos) if v.pos==pos and v.id==vname]
 if len(variants)!=1 or len(variants[0].alleles)!=2 or set(variants[0].alleles)!={count,other}:raise ValueError('reference mismatch '+vname)
 v=variants[0];idx=v.alleles.index(count);g=[]
 for n in names:
  gt=v.samples[n].get('GT')
  if gt is None or len(gt)!=2 or any(a is None for a in gt):raise ValueError('missing genotype '+vname)
  g.append(sum(a==idx for a in gt))
 vectors.append(g);p=float(b['counted_allele_frequency']);beta=float(b['beta_counted_assuming_saige']);se=float(b['se']);pn=float(b['source_p_normal']);ps=float(b['source_p_spa']);assert 0<p<1 and se>0 and 0<pn<=1 and 0<ps<=1
 sign=1 if beta>=0 else -1
 components.append({'variant':vname,'counted_allele':count,'freq':p,'minor_dosage_count':float(b['minor_dosage_count']),'normal_p':pn,'spa_p':ps,'beta_over_se':beta/se,'normal_p_z':sign*NormalDist().inv_cdf(1-pn/2),'spa_p_z':sign*NormalDist().inv_cdf(1-ps/2)})
x=np.array(vectors,float).T;chrom=[r['chr_GRCh37_assumed'] for r in selected];scale=np.array([math.sqrt(2*c['freq']*(1-c['freq'])) for c in components]);result={}
def calc(z,a,R):return float(a@z/math.sqrt(a@R@a))
for group,mat in [('JPT_n104',x[mask]),('EAS_n504',x)]:
 R=np.corrcoef(mat,rowvar=False);R=np.where(np.equal.outer(chrom,chrom),R,0.);assert np.allclose(np.diag(R),1)
 for mode in ('beta_over_se','normal_p_z','spa_p_z'):
  z=np.array([c[mode] for c in components]);val=calc(z,scale,R);loo=[]
  for j in range(14):
   keep=np.arange(14)!=j;loo.append(calc(z[keep],scale[keep],R[np.ix_(keep,keep)]))
  result[group+'_'+mode]={'z':val,'nominal_p':math.erfc(abs(val)/math.sqrt(2)),'leave_one_out_z_min':min(loo),'leave_one_out_z_max':max(loo)}
with (out/'stage6_bbj_insomnia_matched14_z_components.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(components[0]));w.writeheader();w.writerows(components)
receipt={'markers':prior['markers'],'results':result,'inputs_sha256':{name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in ['unvalidated_22_shared_direction_score.csv','stage6_bbj_insomnia_oriented_candidates.csv','stage6_cross_population_matched14.json']},'reference':'1000 Genomes phase 3; cross-chromosome correlations set to zero; JPT n=104 and EAS n=504, overlapping panels, remote VCF not full-file checksum verified','method':'a_i=sqrt(2p_i(1-p_i)) with source frequency, equal marker weights; a.T z/sqrt(a.T R a). Three conditional choices of marginal z; not genotype-score regression.','limits':'Postselected matched 14 subset, BBJ Ncase=389, unverified sample overlap, source allele convention, external LD, no target data. SPA p-derived z is a sensitivity, not a score-level SPA or proof of calibration.'}
(out/'stage6_bbj_insomnia_matched14_z.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(result,indent=2))

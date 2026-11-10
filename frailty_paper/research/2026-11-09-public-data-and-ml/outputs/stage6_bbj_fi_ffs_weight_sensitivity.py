#!/usr/bin/env python3
"""Conditional partial-score summary-Z comparison, not individual prediction."""
import csv,hashlib,json,math
from pathlib import Path
import numpy as np
import pysam
root=Path(__file__).resolve().parent.parent;out=root/'outputs'
def read(name):
 with (out/name).open(newline='') as f:return list(csv.DictReader(f))
score=read('unvalidated_15_variant_score_prototype.csv');bbj={r['variant']:r for r in read('stage6_bbj_insomnia_oriented_candidates.csv') if r['score']=='prototype15' and r['status']=='exact position and pair'};sel=[r for r in score if r['variant'] in bbj];assert len(sel)==13
with (root/'work/1kg_phase3_panel.tsv').open(newline='') as f:panel=list(csv.DictReader(f,delimiter='\t'))
names=sorted(r['sample'] for r in panel if r['super_pop']=='EAS');jpt={r['sample'] for r in panel if r['pop']=='JPT'};jmask=np.array([name in jpt for name in names]);assert len(names)==504 and sum(jmask)==104
handles={};genos=[];items=[]
for r in sel:
 name=r['variant'];b=bbj[name];chrom=r['chr'];position=int(r['position_GRCh37_inferred']);count=r['counted_allele'];other=r['other_allele'];assert (b['counted_allele'],b['other_allele'])==(count,other)
 if chrom not in handles:handles[chrom]=pysam.VariantFile(f'https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz')
 v=[v for v in handles[chrom].fetch(chrom,position-1,position) if v.pos==position and v.id==name]
 if len(v)!=1 or len(v[0].alleles)!=2 or set(v[0].alleles)!={count,other}:raise ValueError('reference mismatch '+name)
 v=v[0];idx=v.alleles.index(count);g=[]
 for n in names:
  gt=v.samples[n].get('GT')
  if gt is None or len(gt)!=2 or any(a is None for a in gt):raise ValueError('missing reference genotype '+name)
  g.append(sum(a==idx for a in gt))
 genos.append(g);p=float(b['counted_allele_frequency']);se=float(b['se']);assert 0<p<1 and se>0
 items.append({'variant':name,'counted_allele':count,'frequency':p,'expected_minor_case_dosage_under_null':2*389*min(p,1-p),'beta_over_se':float(b['beta_counted_assuming_saige'])/se,'source_p_spa':float(b['source_p_spa']),'source_p_normal':float(b['source_p_normal'])})
x=np.array(genos,float).T;chrom=[r['chr'] for r in sel];sets={'matched13':np.ones(13,bool),'expected_case_minor_10_plus':np.array([c['expected_minor_case_dosage_under_null']>=10 for c in items])};assert sets['expected_case_minor_10_plus'].sum()==12
weight_names=['FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight'];result={}
for group,mat in [('JPT_n104',x[jmask]),('EAS_n504',x)]:
 R=np.corrcoef(mat,rowvar=False);R=np.where(np.equal.outer(chrom,chrom),R,0.);assert np.isfinite(R).all()
 for subset,keep in sets.items():
  rr=R[np.ix_(keep,keep)];chosen=[r for r,k in zip(sel,keep) if k];comps=[r for r,k in zip(items,keep) if k];z=np.array([c['beta_over_se'] for c in comps]);p=np.array([c['frequency'] for c in comps]);base=np.sqrt(2*p*(1-p));weights=[np.array([float(r[w]) for r in chosen])*base for w in weight_names]
  for key,a in zip(weight_names,weights):
   value=float(a@z/math.sqrt(a@rr@a)); result[group+'_'+subset+'_'+key]={'n':len(chosen),'z':value,'nominal_p':math.erfc(abs(value)/math.sqrt(2))}
  cov=float(weights[0]@rr@weights[1]/math.sqrt((weights[0]@rr@weights[0])*(weights[1]@rr@weights[1])));result[group+'_'+subset+'_fi_ffs_approx_z_correlation']={'r':cov}
with (out/'stage6_bbj_fi_ffs_weight_components.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(items[0]));w.writeheader();w.writerows(items)
receipt={'inputs_sha256':{n:hashlib.sha256((out/n).read_bytes()).hexdigest() for n in ['unvalidated_15_variant_score_prototype.csv','stage6_bbj_insomnia_oriented_candidates.csv']},'markers':[c['variant'] for c in items],'results':result,'method':'a_i=weight_i sqrt(2p_i(1-p_i)); beta/SE oriented to counted allele; Z=a.T z/sqrt(a.T R a); R within-chromosome reference LD, cross-chromosome zero','limits':'BBJ 389 cases, postselected partial marker sets, low case allele counts; source Allele2 and FFS coding inferred; target LD and cohort overlap uncertain. Comparison cannot establish relative predictive performance, no score-level SPA.'}
(out/'stage6_bbj_fi_ffs_weight_sensitivity.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(result,indent=2))

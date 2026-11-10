#!/usr/bin/env python3
"""Postselected BBJ marker-coverage and low-expected-count sensitivity."""
import csv,hashlib,json,math
from pathlib import Path
import numpy as np
import pysam
root=Path(__file__).resolve().parent.parent;out=root/'outputs'
def read(name):
 with (out/name).open(newline='') as f:return list(csv.DictReader(f))
score=read('unvalidated_22_shared_direction_score.csv');bbj={r['variant']:r for r in read('stage6_bbj_insomnia_oriented_candidates.csv') if r['score']=='shared22' and r['status']=='exact position and pair'}
prior=json.loads((out/'stage6_cross_population_matched14.json').read_text())['markers']
selected=[r for r in score if r['variant'] in bbj and float(bbj[r['variant']]['minor_dosage_count'])>=20]
assert len(selected)==17 and len(bbj)==18
with (root/'work/1kg_phase3_panel.tsv').open(newline='') as f:panel=list(csv.DictReader(f,delimiter='\t'))
names=sorted(r['sample'] for r in panel if r['super_pop']=='EAS');jpt={r['sample'] for r in panel if r['pop']=='JPT'};mask=np.array([n in jpt for n in names]);assert len(names)==504 and sum(mask)==104
handles={};geno=[];items=[]
for r in selected:
 vname=r['variant'];b=bbj[vname];chr=r['chr_GRCh37_assumed'];pos=int(r['position_GRCh37_assumed']);allele=r['counted_allele'];other=r['other_allele'];assert (b['counted_allele'],b['other_allele'])==(allele,other)
 if chr not in handles:handles[chr]=pysam.VariantFile(f'https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz')
 match=[v for v in handles[chr].fetch(chr,pos-1,pos) if v.pos==pos and v.id==vname]
 if len(match)!=1 or len(match[0].alleles)!=2 or set(match[0].alleles)!={allele,other}:raise ValueError('reference mismatch '+vname)
 v=match[0];idx=v.alleles.index(allele);g=[]
 for name in names:
  gt=v.samples[name].get('GT')
  if gt is None or len(gt)!=2 or any(x is None for x in gt):raise ValueError('missing '+vname)
  g.append(sum(a==idx for a in gt))
 geno.append(g);p=float(b['counted_allele_frequency']);z=float(b['beta_counted_assuming_saige'])/float(b['se']);assert 0<p<1
 items.append({'variant':vname,'counted_allele':allele,'freq':p,'minor_frequency':min(p,1-p),'expected_case_minor_dosage_null_approx':2*389*min(p,1-p),'source_beta_over_se':z,'normal_p':float(b['source_p_normal']),'spa_p':float(b['source_p_spa']),'in_matched14':vname in prior})
x=np.array(geno,float).T;chrom=[r['chr_GRCh37_assumed'] for r in selected];sets={'matched14':np.array([c['in_matched14'] for c in items]),'expected_10_or_more_cases':np.array([c['expected_case_minor_dosage_null_approx']>=10 for c in items]),'all_17_after_low_total_count_exclusion':np.ones(17,bool)}
assert sum(sets['matched14'])==14 and sum(sets['expected_10_or_more_cases'])==15
result={}
for group,mat in [('JPT_n104',x[mask]),('EAS_n504',x)]:
 R=np.corrcoef(mat,rowvar=False);R=np.where(np.equal.outer(chrom,chrom),R,0.);assert np.isfinite(R).all()
 for name,keep in sets.items():
  comps=[c for c,k in zip(items,keep) if k];z=np.array([c['source_beta_over_se'] for c in comps]);a=np.array([math.sqrt(2*c['freq']*(1-c['freq'])) for c in comps]);rr=R[np.ix_(keep,keep)];value=float(a@z/math.sqrt(a@rr@a))
  result[group+'_'+name]={'n':len(comps),'z':value,'nominal_p':math.erfc(abs(value)/math.sqrt(2)),'markers':[c['variant'] for c in comps]}
with (out/'stage6_bbj_insomnia_marker_coverage_components.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(items[0]));w.writeheader();w.writerows(items)
receipt={'inputs_sha256':{n:hashlib.sha256((out/n).read_bytes()).hexdigest() for n in ['unvalidated_22_shared_direction_score.csv','stage6_bbj_insomnia_oriented_candidates.csv','stage6_cross_population_matched14.json']},'result':result,'selection':'14 matched earlier sources; 15 with >=10 expected minor copies among 389 cases under null approximation; 17 exclude total minor dosage <20. All postselected; 18th marker excluded due to total minor dosage ~15.','limitations':'Not target score regression; marginal beta/SE, approximate external LD, rare alleles and extreme case-control imbalance, source coding, target LD and overlap uncertain; flags not validated thresholds; no replication or predictive performance.'}
(out/'stage6_bbj_insomnia_marker_coverage_sensitivity.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:{kk:v for kk,v in val.items() if kk!='markers'} for k,val in result.items()},indent=2))

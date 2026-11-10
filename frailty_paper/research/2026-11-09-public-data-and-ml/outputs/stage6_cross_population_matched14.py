#!/usr/bin/env python3
"""Same 14 alleles, two target phenotypes; not a joint-outcome analysis."""
import csv,json,math
from pathlib import Path
import numpy as np
import pysam

root=Path(__file__).resolve().parent.parent;out=root/'outputs'
def read(path,delimiter=','):
    with path.open(newline='') as f:return list(csv.DictReader(f,delimiter=delimiter))
score=read(out/'unvalidated_22_shared_direction_score.csv')
tm={r['variant']:r for r in read(out/'stage6_tpmi_shared22_insomnia_inventory.csv') if r['status']=='browser_lookup_rounded'}
jm={r['variant']:r for r in read(out/'stage6_japanese_knee_candidate_matches.csv') if r['score']=='shared22' and r['exact_rows']=='1'}
selected=[r for r in score if r['variant'] in tm and r['variant'] in jm];assert len(selected)==14 and len(jm)==14
panel=read(root/'work/1kg_phase3_panel.tsv','\t');names=sorted(r['sample'] for r in panel if r['super_pop']=='EAS');pop={r['sample']:r['pop'] for r in panel};assert len(names)==504
is_jpt=np.array([pop[n]=='JPT' for n in names]);assert sum(is_jpt)==104
handles={};genotypes=[];taiwan_z=[];japan_z=[];japan_a=[]
for r in selected:
    name=r['variant'];chrom=r['chr_GRCh37_assumed'];pos=int(r['position_GRCh37_assumed']);count=r['counted_allele'];other=r['other_allele']
    if chrom not in handles:handles[chrom]=pysam.VariantFile(f'https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz')
    matches=[v for v in handles[chrom].fetch(chrom,pos-1,pos) if v.pos==pos and v.id==name]
    if len(matches)!=1 or len(matches[0].alleles)!=2 or set(matches[0].alleles)!={count,other}:raise ValueError('Reference mismatch '+name)
    v=matches[0];idx=v.alleles.index(count);g=[]
    for n in names:
        gt=v.samples[n].get('GT')
        if gt is None or len(gt)!=2 or any(a is None for a in gt):raise ValueError('Missing genotype '+name)
        g.append(sum(a==idx for a in gt))
    genotypes.append(g)
    t=tm[name];assert {t['ref'],t['alt']}=={count,other};tsign=1 if t['alt']==count else -1
    taiwan_z.append(tsign*float(t['insomnia_effect'])/float(t['insomnia_se']))
    j=jm[name];assert {j['source_effect_allele'],j['source_other_allele']}=={count,other};jsign=1 if j['source_effect_allele']==count else -1
    japan_z.append(-jsign*float(j['source_beta'])/float(j['source_standard_error']))
    p=float(j['source_effect_allele_frequency']);japan_a.append(math.sqrt(2*p*(1-p)))
x=np.array(genotypes,dtype=float).T;chrom=[r['chr_GRCh37_assumed'] for r in selected]
def compute(z,a,x):
    R=np.corrcoef(x,rowvar=False);R=np.where(np.equal.outer(chrom,chrom),R,0.);z=np.asarray(z);a=np.asarray(a)
    value=float(a@z/np.sqrt(a@R@a));return {'weighted_z':value,'nominal_two_sided_p':math.erfc(abs(value)/math.sqrt(2))}
taiwan=compute(taiwan_z,x.std(axis=0,ddof=1),x);japan=compute(japan_z,japan_a,x[is_jpt,:])
data={'markers':[r['variant'] for r in selected],'taiwan_insomnia':taiwan,'japan_lower_strength':japan,'limitations':'Same 14 SNPs but different outcomes, ancestry, LD, genotype scaling, sample sizes and study design; postselected, rounded Taiwanese values and conditional ALT direction. Not a joint phenotype, score prediction, replication of full score, genetic correlation, shared causal variant or formal heterogeneity test.'}
(out/'stage6_cross_population_matched14.json').write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(data,indent=2))

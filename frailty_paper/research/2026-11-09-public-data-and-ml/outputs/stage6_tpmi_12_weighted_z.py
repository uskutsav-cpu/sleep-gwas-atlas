#!/usr/bin/env python3
"""Conditional exploratory 12-marker summary screen and reference coverage."""
import csv,hashlib,json,math
from pathlib import Path
import numpy as np
import pysam

root=Path(__file__).resolve().parent.parent;out=root/'outputs'
def read(path,sep=','):
    with path.open(newline='') as f:return list(csv.DictReader(f,delimiter=sep))
score_path=out/'unvalidated_15_variant_score_prototype.csv'; inv_path=out/'stage6_tpmi_insomnia_candidate_inventory.csv'
score=read(score_path);lookup={r['variant']:r for r in read(inv_path)}
names=sorted(r['sample'] for r in read(root/'work/1kg_phase3_panel.tsv','\t') if r['super_pop']=='EAS')
assert len(score)==15 and len(names)==504
handles={};vectors=[]
for r in score:
    chrom=r['chr'];pos=int(r['position_GRCh37_inferred'])
    if chrom not in handles:
        handles[chrom]=pysam.VariantFile(f'https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz')
    m=[v for v in handles[chrom].fetch(chrom,pos-1,pos) if v.pos==pos and v.id==r['variant']]
    if len(m)!=1 or len(m[0].alleles)!=2 or set(m[0].alleles)!={r['counted_allele'],r['other_allele']}:raise ValueError('Unmatched reference '+r['variant'])
    rec=m[0];a=rec.alleles.index(r['counted_allele']);v=[]
    for name in names:
        gt=rec.samples[name].get('GT')
        if gt is None or len(gt)!=2 or any(y is None for y in gt):raise ValueError('Missing reference genotype')
        v.append(sum(y==a for y in gt))
    vectors.append(v)
x=np.array(vectors,dtype=float).T; mask=np.array([lookup[r['variant']]['status']=='browser_lookup_rounded' for r in score]);assert mask.sum()==12
xsub=x[:,mask];sd=xsub.std(axis=0,ddof=1);assert np.all(sd>0)
R=np.corrcoef(xsub,rowvar=False); chroms=[r['chr'] for r,keep in zip(score,mask) if keep];R=np.where(np.equal.outer(chroms,chroms),R,0.)
selected=[r for r,keep in zip(score,mask) if keep]
weights=['FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight'];results={}
for weight in weights:
    w=np.array([float(r[weight]) for r in score]);full=x@w;subset=x[:,mask]@w[mask]
    result={'reference_full_vs_subset_correlation':float(np.corrcoef(full,subset)[0,1]),'reference_full_score_sd':float(full.std(ddof=1)),'reference_subset_score_sd':float(subset.std(ddof=1))}
    for phenotype in ('insomnia','persistent'):
        z=[]
        for r in selected:
            item=lookup[r['variant']]; ref,alt=item['browser_pair'].split('/')
            assert {ref,alt}=={r['counted_allele'],r['other_allele']}
            orient=1 if r['counted_allele']==alt else -1
            z.append(orient*float(item[phenotype+'_effect'])/float(item[phenotype+'_se']))
        z=np.asarray(z);a=w[mask]*sd;stat=float(a@z/np.sqrt(a@R@a));result[phenotype]={'weighted_z':stat,'nominal_two_sided_p':math.erfc(abs(stat)/math.sqrt(2))}
    results[weight]=result
receipt={'method':'a.T z / sqrt(a.T R a), a=score weight times EAS reference dosage SD; cross-chromosome correlations zero','assumption':'browser effect is ALT; conditional for these variants','score_sha256':hashlib.sha256(score_path.read_bytes()).hexdigest(),'inventory_sha256':hashlib.sha256(inv_path.read_bytes()).hexdigest(),'n_markers':12,'n_reference':len(names),'results':results,'limitations':'Postselected subset, rounded browser values and external LD; incomplete and not target score regression or individual validation. FFS weight direction inferred.'}
(out/'stage6_tpmi_12_weighted_z.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(results,indent=2))

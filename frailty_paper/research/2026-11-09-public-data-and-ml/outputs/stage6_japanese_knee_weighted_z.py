#!/usr/bin/env python3
"""Post hoc partial-score marginal-Z screen against knee strength."""
import csv,hashlib,json,math
from pathlib import Path
import numpy as np
import pysam

root=Path(__file__).resolve().parent.parent; out=root/'outputs'
def read(name,sep=','):
    with (out/name).open(newline='') as f:return list(csv.DictReader(f,delimiter=sep))
receipt=json.loads((out/'stage6_japanese_knee_candidate_receipt.json').read_text());assert receipt['checksum_matches']
lookup={(r['score'],r['variant']):r for r in read('stage6_japanese_knee_candidate_matches.csv')}
with (root/'work/1kg_phase3_panel.tsv').open(newline='') as f:
    samples=sorted(r['sample'] for r in csv.DictReader(f,delimiter='\t') if r['pop']=='JPT')
assert len(samples)==104
definitions=[('prototype15','unvalidated_15_variant_score_prototype.csv','chr','position_GRCh37_inferred',['FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight']),('shared22','unvalidated_22_shared_direction_score.csv','chr_GRCh37_assumed','position_GRCh37_assumed',['weight'])]
handles={};results=[];component=[]
for label,filename,ch,pos,weight_names in definitions:
    score=read(filename); genotype=[]
    for r in score:
        chrom=r[ch]; position=int(r[pos])
        if chrom not in handles:handles[chrom]=pysam.VariantFile(f'https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz')
        matches=[v for v in handles[chrom].fetch(chrom,position-1,position) if v.pos==position and v.id==r['variant']]
        if len(matches)!=1 or len(matches[0].alleles)!=2 or set(matches[0].alleles)!={r['counted_allele'],r['other_allele']}:raise ValueError('Reference variant mismatch: '+r['variant'])
        v=matches[0];ai=v.alleles.index(r['counted_allele']);g=[]
        for name in samples:
            gt=v.samples[name].get('GT')
            if gt is None or len(gt)!=2 or any(a is None for a in gt):raise ValueError('Missing reference genotype')
            g.append(sum(a==ai for a in gt))
        genotype.append(g)
    x=np.array(genotype,dtype=float).T; mask=np.array([lookup[(label,r['variant'])]['exact_rows']=='1' for r in score]);assert mask.sum()==({'prototype15':12,'shared22':14}[label])
    subset=[r for r,k in zip(score,mask) if k];xsub=x[:,mask];sd=xsub.std(axis=0,ddof=1);assert np.all(sd>0)
    R=np.corrcoef(xsub,rowvar=False);chroms=[r[ch] for r in subset];R=np.where(np.equal.outer(chroms,chroms),R,0.)
    z=[];scale=[]
    for r in subset:
        d=lookup[(label,r['variant'])];eff=d['source_effect_allele'];freq=float(d['source_effect_allele_frequency']);se=float(d['source_standard_error']);assert se>0 and 0<freq<1
        assert {eff,d['source_other_allele']}=={r['counted_allele'],r['other_allele']}
        orientation=1 if eff==r['counted_allele'] else -1
        # Source outcome = higher strength. Flip to a lower-strength direction.
        z.append(-orientation*float(d['source_beta'])/se);scale.append(math.sqrt(2*freq*(1-freq)))
    z=np.array(z);scale=np.array(scale)
    for weight in weight_names:
        w=np.array([float(r[weight]) for r in score]);a=w[mask]*scale
        stat=float(a@z/np.sqrt(a@R@a));full=x@w;partial=x[:,mask]@w[mask]
        loo=[]
        for j in range(mask.sum()):
            keep=np.arange(mask.sum())!=j;aa=a[keep];val=float(aa@z[keep]/np.sqrt(aa@R[np.ix_(keep,keep)]@aa));loo.append(val)
            component.append({'score':label,'weight':weight,'variant':subset[j]['variant'],'oriented_lower_strength_z':float(z[j]),'leave_one_out_z':val})
        results.append({'score':label,'weight':weight,'n_used':int(mask.sum()),'n_total':len(score),'lower_strength_weighted_z':stat,'nominal_two_sided_p':math.erfc(abs(stat)/math.sqrt(2)),'leave_one_out_min':min(loo),'leave_one_out_max':max(loo),'jpt_full_vs_partial_score_r':float(np.corrcoef(full,partial)[0,1])})
with (out/'stage6_japanese_knee_weighted_z_components.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(component[0]));w.writeheader();w.writerows(component)
data={'source_md5':receipt['observed_md5'],'candidate_table_sha256':hashlib.sha256((out/'stage6_japanese_knee_candidate_matches.csv').read_bytes()).hexdigest(),'reference':'1000 Genomes phase 3 JPT n=104; within-chromosome dosage LD, cross-chromosome zero','method':'a=weight*sqrt(2*target_EAF*(1-target_EAF)); z = negative counted-allele beta / SE; statistic a.T z/sqrt(a.T R a)','results':results,'limitations':'Partial scores, postselected, unharmonized source/strand and multiallelic caveats, FFS direction inferred, small external LD panel, source sample overlap not exhaustively verified; not score regression or individual prediction.'}
(out/'stage6_japanese_knee_weighted_z.json').write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(results,indent=2))

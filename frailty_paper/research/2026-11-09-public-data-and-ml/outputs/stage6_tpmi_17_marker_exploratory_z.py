#!/usr/bin/env python3
"""Conditional, approximate weighted-Z screen; not target score regression."""
import csv, hashlib, json, math
from pathlib import Path
import numpy as np
import pysam

root=Path(__file__).resolve().parent.parent; out=root/'outputs'
def read(path,sep=','):
    with path.open(newline='') as f: return list(csv.DictReader(f,delimiter=sep))
score_path=out/'unvalidated_22_shared_direction_score.csv'
inventory_path=out/'stage6_tpmi_shared22_insomnia_inventory.csv'
score={r['variant']:r for r in read(score_path)}
markers=[r for r in read(inventory_path) if r['status']=='browser_lookup_rounded']
samples=sorted(r['sample'] for r in read(root/'work/1kg_phase3_panel.tsv','\t') if r['super_pop']=='EAS')
assert len(markers)==17 and len(samples)==504
handles={}; vectors=[]
for r in markers:
    s=score[r['variant']]; chrom=s['chr_GRCh37_assumed']; pos=int(s['position_GRCh37_assumed'])
    if chrom not in handles:
        handles[chrom]=pysam.VariantFile(f'https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz')
    records=[v for v in handles[chrom].fetch(chrom,pos-1,pos) if v.pos==pos and v.id==r['variant']]
    if len(records)!=1 or len(records[0].alleles)!=2 or set(records[0].alleles)!={r['counted_allele'],r['other_allele']}:
        raise ValueError(f"Reference mismatch: {r['variant']}")
    rec=records[0]; allele=rec.alleles.index(r['counted_allele']); values=[]
    for name in samples:
        gt=rec.samples[name].get('GT')
        if gt is None or len(gt)!=2 or any(x is None for x in gt): raise ValueError('Missing reference genotype')
        values.append(sum(x==allele for x in gt))
    vectors.append(values)
x=np.asarray(vectors,dtype=float).T
sd=x.std(axis=0,ddof=1); assert np.all(sd>0)
R=np.corrcoef(x,rowvar=False)
chroms=[score[r['variant']]['chr_GRCh37_assumed'] for r in markers]
R=np.where(np.equal.outer(chroms,chroms),R,0.)
assert np.allclose(np.diag(R),1)
results={}; components=[]
for phenotype in ('insomnia','persistent'):
    beta=np.array([float(r[phenotype+'_effect']) for r in markers]); se=np.array([float(r[phenotype+'_se']) for r in markers])
    assert np.all(se>0)
    # Conditional assumption: browser beta is per ALT; orient it to counted allele.
    orientation=np.array([1 if r['counted_allele']==r['alt'] else -1 for r in markers])
    z=beta/se*orientation
    value=float(sd@z / np.sqrt(sd@R@sd)); p=math.erfc(abs(value)/math.sqrt(2))
    leave=[]
    for j in range(len(markers)):
        keep=np.arange(len(markers))!=j; a=sd[keep]; rr=R[np.ix_(keep,keep)]
        leave.append(float(a@z[keep]/np.sqrt(a@rr@a)))
        components.append({'phenotype':phenotype,'variant':markers[j]['variant'],'conditional_counted_z':float(z[j]),'reference_dosage_sd':float(sd[j]),'leave_one_out_z':leave[-1]})
    results[phenotype]={'weighted_z':value,'nominal_two_sided_p':p,'leave_one_out_z_min':min(leave),'leave_one_out_z_max':max(leave)}
with (out/'stage6_tpmi_17_marker_exploratory_components.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(components[0]));w.writeheader();w.writerows(components)
receipt={'method':'a.T z / sqrt(a.T R a), a=EAS reference dosage SD, equal per-allele weights; cross-chromosome correlations set to zero','assumption':'Displayed TPMI beta is ALT effect; conditional, not independently verified for these rows','source_score_sha256':hashlib.sha256(score_path.read_bytes()).hexdigest(),'inventory_sha256':hashlib.sha256(inventory_path.read_bytes()).hexdigest(),'n_markers':17,'n_reference':len(samples),'results':results,'limitations':'Postselected 17-marker subset, rounded manual browser values, external LD, potentially mismatched sample sizes and association models, unresolved multiallelic representations and participant overlap; not full score, target regression, calibration or validation.'}
(out/'stage6_tpmi_17_marker_exploratory_z.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(results,indent=2))

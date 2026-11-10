#!/usr/bin/env python3
"""Compare Japanese source alleles and EAF to 1000 Genomes JPT reference."""
import csv,json
from pathlib import Path
import numpy as np
import pysam

root=Path(__file__).resolve().parent.parent;out=root/'outputs'
with (out/'stage6_japanese_knee_candidate_matches.csv').open(newline='') as f:rows=[r for r in csv.DictReader(f) if r['score']!='paper_control' and r['exact_rows']=='1']
with (root/'work/1kg_phase3_panel.tsv').open(newline='') as f:names=sorted(r['sample'] for r in csv.DictReader(f,delimiter='\t') if r['pop']=='JPT')
assert len(rows)==26 and len(names)==104
handles={}; cache={};results=[]
for r in rows:
    key=(r['chromosome'],int(r['position']),r['variant']);chrom,pos,rs=key
    if key not in cache:
        if chrom not in handles:handles[chrom]=pysam.VariantFile(f'https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz')
        m=[v for v in handles[chrom].fetch(chrom,pos-1,pos) if v.pos==pos and v.id==rs]
        if len(m)!=1 or len(m[0].alleles)!=2:raise ValueError('Reference ambiguity '+rs)
        v=m[0]; dosages={}
        for allele in v.alleles:
            idx=v.alleles.index(allele);a=[]
            for sample in names:
                gt=v.samples[sample].get('GT')
                if gt is None or len(gt)!=2 or any(x is None for x in gt):raise ValueError('Missing genotype '+rs)
                a.append(sum(x==idx for x in gt))
            dosages[allele]=float(np.mean(a)/2)
        cache[key]={'ref':v.alleles[0],'alt':v.alleles[1],'frequencies':dosages}
    ref=cache[key];parts=r['source_variant_id'].split('_');assert len(parts)==4 and parts[0]==chrom and int(parts[1])==pos
    sref,salt=parts[2:];effect=r['source_effect_allele'];assert set([sref,salt])==set([r['source_effect_allele'],r['source_other_allele']])
    if effect not in ref['frequencies']:raise ValueError('Alleles unmatched '+rs)
    target=float(r['source_effect_allele_frequency']);panel=ref['frequencies'][effect]
    results.append({'score':r['score'],'variant':rs,'source_ref':sref,'source_alt':salt,'reference_ref':ref['ref'],'reference_alt':ref['alt'],'exact_reference_order':(sref,salt)==(ref['ref'],ref['alt']),'source_effect_allele':effect,'palindromic':set([sref,salt]) in ({'A','T'},{'C','G'}),'source_effect_frequency':target,'jpt_effect_frequency':panel,'absolute_frequency_difference':abs(target-panel)})
with (out/'stage6_japanese_knee_allele_frequency_audit.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
unique={r['variant']:r for r in results}
summary={'n_rows':len(results),'n_unique_variants':len(unique),'n_reference_order_mismatch':sum(not r['exact_reference_order'] for r in unique.values()),'n_palindromic':sum(r['palindromic'] for r in unique.values()),'max_absolute_frequency_difference':max(r['absolute_frequency_difference'] for r in unique.values()),'n_frequency_difference_gt_0_1':sum(r['absolute_frequency_difference']>.1 for r in unique.values()),'limitations':'JPT panel is 104 individuals, not target sample; source may use differing cohorts and imputation. Frequency agreement and REF/ALT matches do not prove source strand or effect semantics. Reference VCF accessed by region without full-file checksum.'}
(out/'stage6_japanese_knee_allele_frequency_audit.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

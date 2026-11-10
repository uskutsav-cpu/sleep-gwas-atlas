#!/usr/bin/env python3
"""Compare the same fixed 22-allele count in two 1000 Genomes panels."""
import csv, hashlib, json
from pathlib import Path
import numpy as np
import pysam

root=Path(__file__).resolve().parent.parent
out=root/'outputs'
score_file=out/'unvalidated_22_shared_direction_score.csv'
with score_file.open(newline='') as f: score=list(csv.DictReader(f))
with (out/'stage6_tpmi_shared22_insomnia_inventory.csv').open(newline='') as f: browser={r['variant']:r['status']=='browser_lookup_rounded' for r in csv.DictReader(f)}
with (root/'work/1kg_phase3_panel.tsv').open() as f: panel=list(csv.DictReader(f,delimiter='\t'))
groups={group:sorted(r['sample'] for r in panel if r['super_pop']==group) for group in ('EAS','EUR')}
assert len(score)==22 and len(groups['EAS'])==504 and len(groups['EUR'])==503
handles={}; arrays={g:[] for g in groups}; rows=[]
for item in score:
    chrom=item['chr_GRCh37_assumed']; pos=int(item['position_GRCh37_assumed']); name=item['variant']
    if chrom not in handles:
        u=f'https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
        handles[chrom]=pysam.VariantFile(u)
    m=[r for r in handles[chrom].fetch(chrom,pos-1,pos) if r.pos==pos and r.id==name]
    if len(m)!=1 or len(m[0].alleles)!=2 or set(m[0].alleles)!={item['counted_allele'],item['other_allele']}:
        raise ValueError(f'{name}: reference record is missing or ambiguous')
    rec=m[0]; ai=rec.alleles.index(item['counted_allele']); row={'variant':name,'counted_allele':item['counted_allele'],'other_allele':item['other_allele'],'tpmi_browser_lookup':browser[name]}
    for g,names in groups.items():
        d=[]
        for sample in names:
            gt=rec.samples[sample].get('GT')
            if gt is None or len(gt)!=2 or any(x is None for x in gt): raise ValueError(f'{name}: missing or non-diploid {g} genotype')
            d.append(sum(x==ai for x in gt))
        a=np.array(d,dtype=float);arrays[g].append(a); af=float(a.mean()/2)
        row.update({f'{g}_counted_frequency':af,f'{g}_minor_frequency':min(af,1-af),f'{g}_dosage_variance':float(a.var(ddof=1))})
    rows.append(row)
with (out/'stage6_shared22_cross_population_frequencies.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary={}
mask=np.array([browser[r['variant']] for r in score]); assert mask.sum()==17
for g in groups:
    a=np.stack(arrays[g],axis=1); full=a.sum(axis=1);partial=a[:,mask].sum(axis=1)
    summary[g]={'n':len(full),'full_count_mean':float(full.mean()),'full_count_sd':float(full.std(ddof=1)),
        '17_marker_count_mean':float(partial.mean()),'17_marker_count_sd':float(partial.std(ddof=1)),
        'full_vs_17_count_correlation':float(np.corrcoef(full,partial)[0,1]),
        'marker_maf_le_1_percent':sum(r[f'{g}_minor_frequency']<=.01 for r in rows)}
receipt={'score_sha256':hashlib.sha256(score_file.read_bytes()).hexdigest(),'source':'UCSC hg19 1000 Genomes phase 3 indexed VCF and sample panel','panel_sizes':{g:len(v) for g,v in groups.items()},'statistics':summary,'limitations':'Reference samples only; no outcomes; score is unvalidated. VCF was not full-file checksum verified. TPMI labels refer to browser results, not cohort genotype availability.'}
(out/'stage6_shared22_cross_population_reference.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(summary,indent=2))

#!/usr/bin/env python3
"""Count selected alleles in the 1000 Genomes EAS panel, not in TPMI."""
import csv, json
from pathlib import Path
import pysam

root=Path(__file__).resolve().parent.parent
with (root/'work/1kg_phase3_panel.tsv').open() as f:
    samples=sorted(r['sample'] for r in csv.DictReader(f,delimiter='\t') if r['super_pop']=='EAS')
assert len(samples)==504
with (root/'outputs/unvalidated_22_shared_direction_score.csv').open() as f:
    score={r['variant']:r for r in csv.DictReader(f)}
ids=['rs34290943','rs9536347','rs9845387','rs12524473','rs17312861']
handles={}; results=[]
for name in ids:
    r=score[name]; chrom=r['chr_GRCh37_assumed']; pos=int(r['position_GRCh37_assumed'])
    if chrom not in handles:
        url=f'https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
        handles[chrom]=pysam.VariantFile(url)
    matches=[v for v in handles[chrom].fetch(chrom,pos-1,pos) if v.pos==pos and v.id==name]
    row={'variant':name,'chr37':chrom,'pos37':pos,'matches':len(matches),'score_alleles':[r['counted_allele'],r['other_allele']]}
    if len(matches)==1:
        record=matches[0]; row['reference_alleles']=list(record.alleles)
        valid=len(record.alleles)==2 and set(record.alleles)==set(row['score_alleles'])
        row['exact_biallelic_match']=valid
        if valid:
            allele=record.alleles.index(r['counted_allele']); dosages=[]
            for sample in samples:
                gt=record.samples[sample].get('GT')
                if gt is not None and len(gt)==2 and all(x is not None for x in gt): dosages.append(sum(x==allele for x in gt))
            total=sum(dosages); chromosomes=2*len(dosages)
            row.update(called_individuals=len(dosages),counted_copies=total,called_chromosomes=chromosomes,
                       counted_frequency=total/chromosomes if chromosomes else None,
                       minor_frequency=min(total,chromosomes-total)/chromosomes if chromosomes else None)
    results.append(row)
out={'source':'UCSC hg19 1000 Genomes phase 3 VCF, EAS samples','limitation':'Reference-panel frequency is not TPMI frequency, imputation quality, inclusion, or evidence of association. Full source VCF not checksum-verified.','records':results}
(root/'outputs/stage6_tpmi_shared22_missing_eas_frequency.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(results,indent=2))

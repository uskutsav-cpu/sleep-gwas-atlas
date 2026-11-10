#!/usr/bin/env python3
"""Reference-panel frequency check; not target cohort QC."""
from pathlib import Path
import json,pandas as pd,numpy as np,pysam
O=Path('outputs');W=Path('work');p=pd.read_csv(W/'1kg_phase3_panel.tsv',sep='\t');samples=sorted(p.loc[p.super_pop=='EAS','sample']);assert len(samples)==504
score=pd.read_csv(O/'unvalidated_15_variant_score_prototype.csv').set_index('variant'); ids=['rs113851554','rs323509','rs78206187']; handles={};res=[]
for rs in ids:
 r=score.loc[rs];ch=str(r['chr']);pos=int(r.position_GRCh37_inferred)
 if ch not in handles:handles[ch]=pysam.VariantFile(f'https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{ch}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz')
 matches=[v for v in handles[ch].fetch(ch,pos-1,pos) if v.pos==pos and v.id==rs];out={'variant':rs,'chrom':ch,'pos37':pos,'records_at_position_and_rsid':len(matches),'project_allele_pair':[r.counted_allele,r.other_allele]}
 if len(matches)==1:
  v=matches[0];out['reference_alleles']=list(v.alleles);valid=len(v.alleles)==2 and set(v.alleles)=={r.counted_allele,r.other_allele};out['exact_biallelic_match']=valid
  if valid:
   idx=v.alleles.index(r.counted_allele);g=[]
   for name in samples:
    gt=v.samples[name].get('GT')
    if gt is not None and len(gt)==2 and all(x is not None for x in gt):g.append(sum(x==idx for x in gt))
   n=len(g);freq=float(np.mean(g)/2) if n else None;out.update(n_called=n,counted_frequency=freq,minor_frequency=min(freq,1-freq) if freq is not None else None)
 res.append(out)
(O/'stage6_tpmi_missing_eas_frequency.json').write_text(json.dumps({'population':'1000 Genomes phase 3 EAS','source':'UCSC hg19 1000 Genomes phase 3','records':res,'limitation':'Reference frequencies are not TPMI cohort frequencies, INFO or proof of availability. VCF not full-file checksum verified.'},indent=2)+'\n');print(json.dumps(res,indent=2))

#!/usr/bin/env python3
import csv,pysam
from pathlib import Path
root=Path(__file__).resolve().parents[1];work=root/'work';u='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';t=pysam.TabixFile(u,index=str(work/'stage6_1000g_chr7_phase3.tbi'));samples=[x for x in t.header if x.startswith('#CHROM')][0].split('\t')[9:];p={x['sample']:x for x in csv.DictReader((work/'1kg_phase3_panel.tsv').open(),delimiter='\t')};eur=[x for x in samples if x in p and p[x]['super_pop']=='EUR'];idx=[i+1 for i,x in enumerate(eur) if p[x]['pop']=='FIN'];assert len(eur)==503 and len(idx)==99;(work/'stage6_foxp2_fin_reference_indices.txt').write_text('\n'.join(map(str,idx))+'\n');print('FIN',len(idx))

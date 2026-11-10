#!/usr/bin/env python3
"""Join exact-allele molecular and phenotype evidence without treating it as colocalization."""
import csv
from collections import Counter
from pathlib import Path
root=Path(__file__).resolve().parent

def read(name):
    with (root/name).open(newline='') as f:
        reader=csv.DictReader(f);return reader.fieldnames,list(reader)

fields, molecular=read('stage5_variant_gene_tissue_evidence_atlas.csv')
_, phenotypes=read('stage5_variant_multitrait_gwas_matrix.csv')
_, queue=read('stage5_score_queue.csv')
idx={r['variant_GRCh38_ref_alt']:r for r in phenotypes}
assert len(idx)==len(phenotypes)==len(queue)
assert set(idx)=={r['variant_key'] for r in queue}
assert {r['variant_GRCh38_ref_alt'] for r in molecular}==set(idx)
traits=['Atkins_FI','Jansen_insomnia','Mak_FI','Mak_FP','Ye_FFS','FinnGen_R13_insomnia']
extra=[f'{t}_{c}' for t in traits for c in ('match_status','beta_alt','se','p')]
assert not set(fields)&set(extra)
out=root/'stage5_variant_gene_tissue_phenotype_integrated.csv'
with out.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields+extra);w.writeheader()
    for row in molecular:
        g=idx[row['variant_GRCh38_ref_alt']]
        row.update({k:g[k] for k in extra});w.writerow(row)
print('rows',len(molecular),'variants',len(idx),'molecular rows with gene',sum(bool(r['gene_id'] or r['gene_symbol']) for r in molecular))
for t in traits:
    print(t,'variants with p value',sum(bool(r[t+'_p']) for r in phenotypes))

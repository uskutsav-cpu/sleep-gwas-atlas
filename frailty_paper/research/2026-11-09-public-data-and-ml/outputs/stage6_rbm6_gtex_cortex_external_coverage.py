#!/usr/bin/env python3
"""Audit published GTEx RBM6 credible-set variants against the coloc input."""
import csv, gzip, hashlib, json
from pathlib import Path
import pandas as pd
root=Path(__file__).resolve().parent.parent
src=root/'work/QTD000171.credible_sets.tsv.gz'
with gzip.open(src,'rt') as f:
    rows=[r for r in csv.DictReader(f,delimiter='\t') if r['gene_id']=='ENSG00000004534']
x=pd.read_csv(root/'outputs/stage6_rbm6_gtex_cortex_aligned.csv')
y=pd.read_csv(root/'outputs/stage6_rbm6_brain_gnomad_aligned.csv')
keys=set(x.variant)
for r in rows: r['in_coloc_input']=str(r['variant'] in keys).lower()
unique={}
for r in rows:
    key=(r['cs_id'],r['variant'])
    if key in unique:
        for field in ('pip','pvalue','beta','se','z','cs_size','in_coloc_input'):
            assert unique[key][field] == r[field], (key,field)
    else:
        unique[key]=r
out=root/'outputs/stage6_rbm6_gtex_cortex_published_cs_coverage.tsv'
with out.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
z=x.merge(y,on=['pos','pair','anchor_effect_allele','anchor_other_allele'],validate='one_to_one')
cs=[]
for cid in sorted({r['cs_id'] for r in rows}):
    v=[r for r in unique.values() if r['cs_id']==cid]; top=max(v,key=lambda r: float(r['pip']))
    assert len(v)==int(v[0]['cs_size'])
    cs.append({'cs_id':cid,'n_variants':len(v),'n_matched':sum(r['variant'] in keys for r in v),'sum_reported_pip':sum(float(r['pip']) for r in v),'sum_reported_pip_matched':sum(float(r['pip']) for r in v if r['variant'] in keys),'top_variant':top['variant'],'top_rsid':top['rsid'],'top_pip':float(top['pip']),'top_in_input':top['variant'] in keys})
receipt={'source_url':'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000015/QTD000171/QTD000171.credible_sets.tsv.gz','source_sha256_local':hashlib.sha256(src.read_bytes()).hexdigest(),'published_checksum_verified':False,'source_rows':len(rows),'sets':cs,'gtex_singlebrain_matched_variants':len(z),'correlation_signed_z_descriptive':float(z.qtl_z_anchor.corr(z.brain_z_anchor)),'note':'Deduplicated identical variant aliases before PIP sums; PIPs describe published credible-set variants only. Correlated variants are not independent observations.'}
assert len(x)==947 and len(z)==501 and len(rows)==52
(root/'outputs/stage6_rbm6_gtex_cortex_external_coverage_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))

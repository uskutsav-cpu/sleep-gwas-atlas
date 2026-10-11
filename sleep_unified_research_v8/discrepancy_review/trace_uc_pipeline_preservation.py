#!/usr/bin/env python3
"""Bounded UC preservation check: first20 harmonized records, before values inspected.
This is NOT a disease-orientation validator; it checks whether source Effect
and Allele1/2 were preserved through the pipeline on sampled retained rows.
"""
import csv,gzip,itertools,json
from pathlib import Path
OUT=Path(__file__).resolve().parent
DATA=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data')
def read(p):
 with gzip.open(p,'rt') as h:yield from csv.DictReader(h,delimiter='\t')
cs=list(itertools.islice(read(DATA/'harmonized/uc.harmonized.tsv.gz'),20))
coords={(c['CHR'],c['BP']) for c in cs};raw={}
for x in read(DATA/'raw/uc.txt.gz'):
 k=tuple(x['MarkerName'].split('_')[0].split(':'))
 if k in coords:raw[k]=x
 if len(raw)==len(coords):break
ids={c['SNP'] for c in cs};ms={x['SNP']:x for x in read(DATA/'munged/uc.sumstats.gz') if x['SNP'] in ids}
rs=[]
for h in cs:
 x=raw[h['CHR'],h['BP']];m=ms.get(h['SNP'],{})
 a=(x['Allele1'].upper(),x['Allele2'].upper());b=(h['A1'],h['A2']);c=(m.get('A1'),m.get('A2'))
 expected_sign=(1 if float(h['BETA'])>0 else -1)* (1 if b==c else -1 if b==c[::-1] else 0)
 rs.append(dict(SNP=h['SNP'],CHR=h['CHR'],BP=h['BP'],raw_A1=a[0],raw_A2=a[1],raw_Effect=x['Effect'],harmonized_A1=b[0],harmonized_A2=b[1],harmonized_BETA=h['BETA'],munged_A1=c[0],munged_A2=c[1],munged_Z=m.get('Z',''),raw_harmonized_identity_beta_match=a==b and float(x['Effect'])==float(h['BETA']),harmonized_munged_oriented_sign_match=bool(m) and expected_sign*float(m['Z'])>0))
with (OUT/'private_uc_pipeline_preservation20.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
s={'sample_rule':'first20 frozen harmonized rows before values inspected','rows':20,'raw_harmonized_identity_beta_match':sum(x['raw_harmonized_identity_beta_match'] for x in rs),'munged_present':sum(bool(x['munged_Z']) for x in rs),'harmonized_munged_oriented_sign_match':sum(x['harmonized_munged_oriented_sign_match'] for x in rs),'certifies_disease_direction':False}
(OUT/'uc_pipeline_preservation_summary.json').write_text(json.dumps(s,indent=2)+'\n');print(s)

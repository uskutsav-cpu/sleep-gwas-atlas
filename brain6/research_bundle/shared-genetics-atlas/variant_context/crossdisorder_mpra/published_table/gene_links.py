#!/usr/bin/env python3
import csv,hashlib,json,zipfile
from pathlib import Path
import openpyxl
p=Path(__file__).resolve().parent
xlsx=p/'table_s4.xlsx'
with zipfile.ZipFile(xlsx) as z:assert z.testzip() is None
with (p/'candidate_comparison.tsv').open() as f: candidates={x['rsid']:x for x in csv.DictReader(f,delimiter='\t') if x['allelic_class']=='emVar'}
assert set(candidates)=={'rs325485','rs6874138'}
w=openpyxl.load_workbook(xlsx,read_only=True,data_only=True)
assert w.sheetnames==['Description','Table S4A','Table S4B','Table S4C','Table S4D','Table S4E']
targets={'rs325485':{'pos37':103995368,'pos38':104659667,'alleles':{'A','G'}},'rs6874138':{'pos37':103899596,'pos38':104563895,'alleles':{'T','C'}}}
results=[];counts={};mismatches=[]
for sheet in ('Table S4A','Table S4B','Table S4C'):
 it=w[sheet].iter_rows(values_only=True);headers=next(it);rows=0
 for row in it:
  if not any(x is not None for x in row):continue
  rows+=1;r=dict(zip(headers,row));rs=r['RSID']
  if rs not in targets:continue
  target=targets[rs];chrom=str(r['CHR']).replace('chr','');pos=int(r['BP']);build='pos37' if sheet=='Table S4C' else 'pos38'
  valid=chrom=='5' and pos==target[build] and (sheet=='Table S4C' or {r['A1'],r['A2']}==target['alleles'])
  item={'sheet':sheet,'rsid':rs,'chromosome':chrom,'position':pos,'genome_build':'GRCh37' if sheet=='Table S4C' else 'GRCh38','matched':valid}
  for key in ('ENSG','Method','eQTL_HGNCName','eQTL_GeneInfo','eQTL_P','eQTL_FDR','iSNP','mDisorder'): 
   if key in r:item[key]=r[key]
  (results if valid else mismatches).append(item)
 counts[sheet]=rows
receipt={'source_url':'https://ars.els-cdn.com/content/image/1-s2.0-S0092867424014351-mmc5.xlsx','source_sha256':hashlib.sha256(xlsx.read_bytes()).hexdigest(),'source_bytes':xlsx.stat().st_size,'xlsx_zip_integrity':'passed','worksheet_rows':counts,'target_variants':sorted(targets),'matched_rows':results,'identifier_matches_with_coordinate_or_allele_mismatch':mismatches,'limitations':'S4C has no allele columns. Results are selected author mappings, not exhaustive molecular tests, independent sleep–ADHD replication, or causal gene mediation. QTL and MPRA effect directions remain unverified.'}
(p/'gene_link_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

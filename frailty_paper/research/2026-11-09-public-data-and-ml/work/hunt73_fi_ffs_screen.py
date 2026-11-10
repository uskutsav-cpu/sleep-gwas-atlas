#!/usr/bin/env python3
"""Compare targeted HUNT SNP identifiers and allele labels to cached FI and FFS summaries."""
import csv,gzip,collections,hashlib,json
from pathlib import Path
import openpyxl
hunt=collections.defaultdict(list)
for r in openpyxl.load_workbook('work/hunt_circadian_supp.xlsx',read_only=True,data_only=True).active.iter_rows(min_row=2,values_only=True):
 if isinstance(r[0],str) and r[0].startswith('rs'): hunt[r[0]].append(r)
assert len(hunt)==73
base=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1')
sources={'FI_HM3':base/'analysis-workspace/frailty_v1/harmonized/frailty.harmonized.tsv.gz','FFS':base/'frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv'}
rows=[];receipt={}
for name,path in sources.items():
 found=collections.defaultdict(list);n=0;openfile=gzip.open if path.suffix=='.gz' else open
 with openfile(path,'rt',newline='') as f:
  for r in csv.DictReader(f,delimiter='\t'):
   n+=1
   if r['SNP'] in hunt:found[r['SNP']].append(r)
 h=hashlib.sha256()
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
 for rs,variants in sorted(hunt.items()):
  a=variants[0]; pairs=[x for x in found.get(rs,[]) if {x['A1'].upper(),x['A2'].upper()}=={str(a[3]).upper(),str(a[4]).upper()}]
  for x in pairs or [{}]:
   rows.append(dict(source=name,rsid=rs,hunt_gene=a[2],hunt_allele_pair=f'{a[3]}/{a[4]}',source_rs_rows=len(found.get(rs,[])),pair_rows=len(pairs),source_chr=x.get('CHR',x.get('CHROM','')),source_pos=x.get('BP',x.get('POS','')),source_a1=x.get('A1',''),source_a2=x.get('A2',''),source_beta=x.get('BETA',''),source_p=x.get('P','')))
 receipt[name]={'rows_scanned':n,'sha256':h.hexdigest(),'unique_ids_found':len(found),'matched_pairs':sum(bool(r['pair_rows']) for r in rows if r['source']==name)}
 print(name,json.dumps(receipt[name]),flush=True)
with open('outputs/stage6_hunt73_fi_ffs_allele_screen.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
with open('outputs/stage6_hunt73_fi_ffs_receipt.json','w') as f:json.dump(receipt,f,indent=2)

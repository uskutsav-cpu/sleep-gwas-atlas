#!/usr/bin/env python3
"""Post-selected exact-ID and coordinate lookup, retaining all candidate records."""
import csv, gzip, hashlib, json
from pathlib import Path

root = Path(__file__).resolve().parents[2]
here = Path(__file__).resolve().parent
info = json.loads((root / 'source_provenance/raw_source_hashes.json').read_text())['FinnGen_R13_insomnia_GRCh38']
sources = {
 'FI': (Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/data/raw/frailty.txt.gz'), 'chromosome', 'base_pair_location', {94666601,94670582}),
 'FFS': (Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv'), 'CHROM', 'POS', {94666601,94670582}),
 'R13_insomnia': (Path(info['path']), '#chrom', 'pos', {94272825}),
}
results={}
for label,(p,chr_col,pos_col,positions) in sources.items():
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
 digest=h.hexdigest()
 if label=='R13_insomnia':assert digest==info['sha256']
 opener=gzip.open if p.suffix=='.gz' else open
 hits=[];n=0
 with opener(p,'rt',newline='') as f:
  rd=csv.DictReader(f,delimiter='\t')
  for row in rd:
   n+=1
   if row[chr_col]=='12' and int(row[pos_col]) in positions:hits.append(row)
 results[label]={'source':str(p),'sha256':digest,'rows_scanned':n,'coordinate_rows':hits}
 print(label, n, len(hits), flush=True)
(here/'receipt.json').write_text(json.dumps(results,indent=2)+'\n')

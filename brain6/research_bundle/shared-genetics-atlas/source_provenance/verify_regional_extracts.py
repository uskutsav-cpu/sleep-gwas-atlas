#!/usr/bin/env python3
"""Re-extract the saved windows from the currently observed source files, without altering them."""
import csv,gzip,hashlib,io,json,pathlib,sys
base=pathlib.Path(__file__).resolve().parent.parent
kind=sys.argv[1]
files=[base/'frailty_coloc'/f'{kind}_extraction.json',base/'frailty_coloc_chr5'/f'{kind}_extraction.json']
meta=[json.loads(p.read_text()) for p in files]
assert meta[0]['path']==meta[1]['path']
keys=('CHR','BP') if kind=='FI' else ('#chrom','pos')
rows=[[],[]]; count=0
with gzip.open(meta[0]['path'],'rt') as f:
 for row in csv.DictReader(f,delimiter='\t'):
  count+=1; chrom=row[keys[0]]
  for i,m in enumerate(meta):
   if chrom in m['regions'] and m['regions'][chrom][0]<=int(row[keys[1]])<=m['regions'][chrom][1]: rows[i].append(row)
results=[]
for m,part,p in zip(meta,rows,files):
 s=io.StringIO(newline='')
 if part:
  w=csv.DictWriter(s,fieldnames=part[0],delimiter='\t');w.writeheader();w.writerows(part)
 digest=hashlib.sha256(s.getvalue().encode()).hexdigest()
 results.append({'extraction':str(p.relative_to(base)),'source_rows':count,'region_rows':len(part),'recomputed_extract_sha256':digest,'saved_extract_sha256':m['extraction_sha256'],'matches_saved':digest==m['extraction_sha256'],'counts_match':count==m['n_source_rows'] and len(part)==m['n_extracted']})
 print(kind,results[-1],flush=True)
(base/'source_provenance'/f'{kind}_reextraction_check.json').write_text(json.dumps(results,indent=2)+'\n')
if not all(r['matches_saved'] and r['counts_match'] for r in results):sys.exit(1)

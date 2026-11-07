import csv, json, hashlib
from pathlib import Path
import openpyxl
root=Path('/Volumes/Extreme SSD/brain6-matoba-mpra')
out=Path(__file__).resolve().parent
ids=['rs62264764','rs62264767','rs2431108','rs77960','rs416223','rs1476535','rs4916723']
# Candidate coordinate probes, using previously verified genome-assembly annotations.
coords={('chr3',117920728),('chr3',117923158),('chr3',117639575),('chr3',117642005),('chr5',104612267),('chr5',104628884),('chr5',104655775),('chr5',88558577),('chr7',114430980)}
regions=[('chr3',117100000,118500000),('chr5',103400000,105100000),('chr5',87000000,89000000),('chr7',113500000,115000000)]
results=[]; meta=[]; regional=[]
for n in [2,3]:
 p=root/f'table{n}.xlsx'; w=openpyxl.load_workbook(p,read_only=True,data_only=True)
 for s in w:
  if n==2 and s.title not in ['a. Regulatory activity','b. allelic effects']: continue
  if n==3 and s.title=='Legends': continue
  count=0; matches=[]; regional_ids={r:set() for r in regions}
  for row in s.iter_rows(values_only=True):
   if not any(x is not None for x in row): continue
   count+=1
   for i,x in enumerate(row):
    if str(x) in ids: matches.append((str(x),'rsID',str(row)))
   if n==2 and count>1:
    ci,pi=(8,9) if s.title.startswith('a.') else (10,11)
    if len(row)>pi and (row[ci],row[pi]) in coords: matches.append((f'{row[ci]}:{row[pi]}','coordinate',str(row)))
   if count>1 and (n==2 or s.title.startswith(('a. Regulatory','c. Regulatory'))):
    ci,pi=((8,9) if s.title.startswith('a.') else (10,11)) if n==2 else (8,9)
    chrom=str(row[ci]); chrom=chrom if chrom.startswith('chr') else 'chr'+chrom
    try: pos=int(row[pi])
    except (TypeError,ValueError): continue
    for reg in regions:
     if chrom==reg[0] and reg[1]<=pos<=reg[2]: regional_ids[reg].add(str(row[0] if n==2 else row[1]))
  if n==2 or s.title.startswith(('a. Regulatory','c. Regulatory')):
   for reg,found in regional_ids.items(): regional.append({'table':n,'sheet':s.title,'region':list(reg),'unique_ids':len(found),'ids':sorted(found)})
  results.append({'table':n,'sheet':s.title,'nonempty_rows_including_header':count,'matches':len(matches),'details':matches})
 meta.append({'table':n,'filename':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(out/'lookup_receipt.json').write_text(json.dumps({'sources':meta,'candidate_rsids':ids,'coordinate_probes':sorted([list(x) for x in coords]),'results':results,'regional_coverage':regional},indent=2)+'\n')
print(json.dumps(results,indent=2))

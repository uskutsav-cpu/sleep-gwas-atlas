import csv, hashlib, io, json, re, zipfile
from pathlib import Path
import requests, openpyxl

out=Path(__file__).resolve().parent
url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10873509/supplementaryFiles'
r=requests.get(url,timeout=120);r.raise_for_status();z=zipfile.ZipFile(io.BytesIO(r.content))
receipt={'url':url,'archive_bytes':len(r.content),'archive_sha256':hashlib.sha256(r.content).hexdigest(),'source_article':'https://doi.org/10.1038/s41598-024-54302-2','assembly':'not independently resolved for reporter coordinates','members':{},'queries':{'rs62264764_grch38_1based':117920728,'rs62264767_grch38_1based':117923158,'alphagenome_bed_window_grch38':[117920477,117920978],'rs62264764_grch37_1based':117639575,'rs62264767_grch37_1based':117642005,'grch37_wide_proximity_bed':[117639074,117640075]}}
rowsout=[]
for no in (2,3):
 name=f'41598_2024_54302_MOESM{no}_ESM.xlsx';info=z.getinfo(name);raw=z.read(info)
 w=openpyxl.load_workbook(io.BytesIO(raw),read_only=True,data_only=True);sheet=w.worksheets[1];iterator=sheet.iter_rows(values_only=True);header=list(next(iterator));n=0;valid=0;chr3=0;nearest=[];exact=[];overlap=[]
 for row in iterator:
  n+=1
  if no==2:
   m=re.fullmatch(r'(?:chr)?(\d+):(\d+)-(\d+)',str(row[0]))
   if not m:continue
   chrom,start,end=m.group(1),int(m.group(2)),int(m.group(3))
  else:
   try:chrom=str(row[0]).removeprefix('chr');start=int(row[1]);end=int(row[2])
   except (ValueError,TypeError):continue
  valid+=1
  if chrom!='3':continue
  chr3+=1
  flags=[]
  for nameq,pos in [('rs62264764_grch38',117920728),('rs62264767_grch38',117923158),('rs62264764_grch37',117639575),('rs62264767_grch37',117642005)]:
   if start<=pos-1<end:flags.append(nameq)
  if start<117920978 and end>117920477:flags.append('window_grch38')
  if start<117640075 and end>117639074:flags.append('wide_proximity_grch37')
  distance=max(start-117920727,117920727-end+1,0)
  nearest.append((distance,start,end))
  if flags:
   rec={'workbook':no,'row':n+1,'flags':','.join(flags),'source_row':list(row)};rowsout.append(rec);overlap.append(rec)
 receipt['members'][str(no)]={'name':info.filename,'bytes':len(raw),'crc32':f'{info.CRC:08x}','sha256':hashlib.sha256(raw).hexdigest(),'sheet':sheet.title,'header':header,'data_rows':n,'parseable_coordinate_rows':valid,'chr3_rows':chr3,'matched_rows':len(overlap),'nearest_to_candidate_assuming_same_build_and_BED':sorted(nearest)[:3]}
(out/'receipt.json').write_text(json.dumps(receipt,indent=2,default=str)+'\n')
(out/'matched_rows.json').write_text(json.dumps(rowsout,indent=2,default=str)+'\n')
print(json.dumps(receipt,indent=2))

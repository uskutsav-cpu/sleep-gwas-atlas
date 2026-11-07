import csv,hashlib,json
from pathlib import Path
import openpyxl
root=Path(__file__).resolve().parent
src=Path('/Volumes/Extreme SSD/brain6-neuropsychiatric-mpra/data5.xlsx')
w=openpyxl.load_workbook(src,read_only=True,data_only=True)
target='chr5_103981726';records=[];raw=[]
for sheet in w:
 if sheet.title in ('P-NPC','correlation_samples'):continue
 if sheet.title=='raw_counts':
  for row in sheet.iter_rows(values_only=True):
   if row[1] and str(row[1]).startswith('chr5.103981726.'):raw.append(str(row[1]))
  continue
 for row in sheet.iter_rows(values_only=True):
  if len(row)>1 and row[1]==target:
   records.append(dict(cell=sheet.title,statistic=row[2],p=row[3],fdr=row[4],logFC=row[5]))
  elif row[0]==target:
   records.append(dict(cell=sheet.title,statistic=row[1],p=row[2],fdr=row[3],logFC=row[8]))
assert len(records)==13 and len(raw)==20 and len({r['cell'] for r in records})==13
with (root/'cell_results.tsv').open('w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=['cell','statistic','p','fdr','logFC'],delimiter='\t');wr.writeheader();wr.writerows(records)
(root/'receipt.json').write_text(json.dumps({'source_file':str(src),'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'variant_coordinate':'GRCh37 chr5:103981726','coordinate_match':'rs40465; independent existing atlas GRCh38 5:104646025 T>G','raw_construct_ids':raw,'excluded_sheet':'P-NPC lacks row identifiers in the published worksheet','matched_cell_count':len(records),'min_reported_fdr':min(r['fdr'] for r in records)},indent=2)+'\n')
print('matched cells',len(records),'minimum reported FDR',min(r['fdr'] for r in records))

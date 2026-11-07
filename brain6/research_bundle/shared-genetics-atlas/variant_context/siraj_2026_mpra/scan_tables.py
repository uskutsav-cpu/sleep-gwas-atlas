import csv,zipfile,gzip,json,hashlib
from pathlib import Path
import openpyxl
p=Path(__file__).resolve().parent;ssd=Path('/Volumes/Extreme SSD/brain6-large-mpra');src=p.parent/'rs40465_mpra/ranked_screen.tsv'
with src.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
wanted={r['grch37']+':'+r['grch38'].split(':')[-1].replace('>',':'):r for r in rows if r.get('grch37')}
assert len(wanted)==425
found=[];n=0
with zipfile.ZipFile(ssd/'tables.zip') as z:
 with gzip.open(z.open('tables/stable2.txt.gz'),'rt') as f:
  for r in csv.DictReader(f,delimiter='\t'):
   n+=1
   if r['variant'] in wanted:found.append({'rsid':wanted[r['variant']]['rsid'],'region':wanted[r['variant']]['region'],'variant37':r['variant'],'duplicate':r['duplicate']})
with (p/'design_matches.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['rsid','region','variant37','duplicate'],delimiter='\t');w.writeheader();w.writerows(found)
w=openpyxl.load_workbook(ssd/'table7.xlsx',read_only=True,data_only=True);sheet=w.active;results=[];total=0
for row in sheet.iter_rows(values_only=True):
 total+=1
 if row[0] in wanted:results.append([wanted[row[0]]['rsid']]+list(row))
with (p/'cell_specific_matches.tsv').open('w',newline='') as f:
 writer=csv.writer(f,delimiter='\t');writer.writerow(['rsid','variant','cell_type','cohort','posterior_mean_activity','sd_activity','lfsr_activity','number_celltypes_active','cell_type_specific_active','posterior_mean_skew','sd_skew','lfsr_skew','number_celltypes_skew','cell_type_specific_skew']);writer.writerows(results)
sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest()
receipt={'design_rows':n,'matched_rows':len(found),'cell_table_rows_including_header':total,'cell_table_matches':len(results),'input_sha256':{str(x):sha(x) for x in [src,ssd/'tables.zip',ssd/'table7.xlsx']},'output_sha256':{x.name:sha(x) for x in [p/'design_matches.tsv',p/'cell_specific_matches.tsv']},'build':'candidate coordinates converted from GRCh38 to GRCh37 using previously documented chain; verify source build before effect interpretation'}
(p/'table_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

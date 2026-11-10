import csv,io,zipfile,hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1];p=root/'work/tenk10k_corrected_caqtl_eqtl.zip';hits=[];n=0
with zipfile.ZipFile(p) as z:
 for name in z.namelist():
  if not name.endswith('/chr5.csv'):continue
  for row in csv.DictReader(io.TextIOWrapper(z.open(name))):
   n+=1
   if row.get('gene')=='ENSG00000113658':row['source_file']=name;hits.append(row)
assert hits
with (root/'outputs/stage6_tenk10k_smad5_molecular_coloc.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(hits[0]));w.writeheader();w.writerows(hits)
print('sha256',hashlib.sha256(p.read_bytes()).hexdigest(),'chr5 rows',n,'SMAD5 rows',len(hits),'H4 >=0.8',sum(float(r['PP.H4.abf'])>=0.8 for r in hits))

import csv,gzip,hashlib,json,math
from pathlib import Path
import requests
p=Path(__file__).resolve().parent;dest=Path('/Volumes/Brain6NHANES/neale-falls-20261028/falls.ldsc.bgz');assert not dest.exists()
u='https://www.dropbox.com/s/2irydbbyi4yi6bb/2296.ldsc.imputed_v3.both_sexes.tsv.bgz?dl=1';r=requests.get(u,timeout=(30,90),stream=True);r.raise_for_status();h=hashlib.sha256();n=0
with dest.open('xb') as f:
 for c in r.iter_content(1024*1024):
  if c:n+=len(c);assert n<9000000;h.update(c);f.write(c)
ids=set();dups=0;invalid=0;N=set();rows=0
with gzip.open(dest,'rt') as f:
 x=csv.DictReader(f,delimiter='\t');fields=x.fieldnames
 for q in x:
  rows+=1;dups+=q['SNP'] in ids;ids.add(q['SNP']);z=float(q['Z']);nval=float(q['N']);N.add(nval)
  if not math.isfinite(z) or not 0<nval<1e9 or q['A1'] not in 'ACGT' or q['A2'] not in 'ACGT':invalid+=1
receipt={'original_url':u,'publication_page':'https://nealelab.github.io/UKBB_ldsc/h2_summary_2296.html','path':str(dest),'bytes':n,'sha256':h.hexdigest(),'rows':rows,'duplicate_ids':dups,'invalid_rows':invalid,'sample_sizes':sorted(N),'header':fields,'gzip_read_completed':True,'publisher_checksum_available':False}
(p/'download_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

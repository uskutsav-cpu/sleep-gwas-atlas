import csv,gzip,hashlib,io,json,subprocess,pathlib
url='https://zenodo.org/records/17069328/files/mv-Sarcopenia.gwas-ssf.tsv.gz?download=1'
expected='d9b5c1e7513e209bc57326264c95c356'; variants=set()
for p in ['outputs/unvalidated_15_variant_score_prototype.csv','outputs/expanded_22_candidate_risk_gene_table.csv']:
 with open(p) as f:variants.update(r['variant'] for r in csv.DictReader(f))
class Hashed:
 def __init__(self,r):self.r=r;self.md5=hashlib.md5();self.n=0
 def read(self,n=-1):
  b=self.r.read(n);self.md5.update(b);self.n+=len(b);return b
 def close(self):self.r.close()
 def seekable(self):return False
 def readable(self):return True
proc=subprocess.Popen(['curl','-L','--fail','--max-time','900','-sS',url],stdout=subprocess.PIPE)
with proc.stdout as raw:
 h=Hashed(raw); z=gzip.GzipFile(fileobj=h); text=io.TextIOWrapper(z,encoding='utf8',newline=''); reader=csv.DictReader(text,delimiter='\t')
 rows=[]; n=0
 for row in reader:
  n+=1
  if row['rsid'] in variants:rows.append(row)
 digest=h.md5.hexdigest(); length=h.n
 if proc.wait()!=0:raise ValueError('download failed; no result output')
 if digest!=expected:raise ValueError(f'compressed checksum mismatch {digest}; no result output')
out=pathlib.Path('outputs/stage6_zenodo_sarcopenia_candidates.csv')
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=reader.fieldnames);w.writeheader();w.writerows(rows)
receipt={'url':url,'md5':digest,'bytes':length,'rows':n,'candidates_total':len(variants),'matched_rows':len(rows),'matched_unique':len({r['rsid'] for r in rows}),'columns':reader.fieldnames,'warning':'Source UK Biobank overlaps discovery; metadata, sample size, genome build and latent direction need independent verification. No individual score validation.'}
pathlib.Path('outputs/stage6_zenodo_sarcopenia_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))

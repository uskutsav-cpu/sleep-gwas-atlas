import csv,gzip,collections,json,hashlib
from pathlib import Path
p=Path(__file__).resolve().parent;parent=p.parent;receipt=json.loads((parent/'retrieval_receipt.json').read_text());source=Path(next(x['file'] for x in receipt['members'] if x['cell']=='progenitor'));peak='chr5_104610871_104612270'
def sha(x):
 h=hashlib.sha256()
 with x.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
assert sha(source)==next(x['inner_gzip_sha256'] for x in receipt['members'] if x['cell']=='progenitor')
rows=[];n=0
with gzip.open(source,'rt') as f:
 header=f.readline().split();assert header==['snp','beta_coef','pvalue','peak','EffectAllele']
 for line in f:
  n+=1;x=line.split()
  if x[3]==peak:rows.append(x)
def key(v):
 chrom,pos,a,b=v.split(':');
 if chrom!='chr5' or a not in 'ACGT' or b not in 'ACGT' or a==b:return None
 return (int(pos),tuple(sorted((a,b))))
valid=[key(x[0]) for x in rows];counts=collections.Counter(x for x in valid if x is not None);unique={k for k,v in counts.items() if v==1};out={}
for trait in ('insomnia','adhd'):
 f=p.parents[1]/'chr5_dense_qtl'/f'Testis_ENSG00000251574_{trait}_matched.tsv'
 with f.open() as h:r=list(csv.DictReader(h,delimiter='\t'))
 g=collections.Counter((int(x['pos38']),tuple(sorted((x['ref'],x['alt'])))) for x in r);guniq={k for k,v in g.items() if v==1};out[trait]={'matched_gwas_qtl_rows':len(r),'unique_reference_keys':len(guniq),'overlap_unique_variant_keys':len(unique & guniq),'input_sha256':sha(f)}
result={'peak':peak,'source_header':header,'source_rows_scanned':n,'peak_rows':len(rows),'eligible_unique_keys':len(unique),'duplicate_variant_keys':sum(v>1 for v in counts.values()),'invalid_variant_rows':sum(x is None for x in valid),'source_reported_se':False,'source_reported_frequency':False,'gwascross':out,'source_sha256':sha(source),'note':'Overlap counts only; selected Model1A support does not establish complete cis coverage or statistical validity.'}
(p/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

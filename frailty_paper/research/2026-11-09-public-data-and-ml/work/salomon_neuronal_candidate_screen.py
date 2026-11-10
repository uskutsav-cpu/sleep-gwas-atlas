import csv,gzip,json,hashlib,collections
from pathlib import Path
root=Path('outputs'); f=Path('work/salomon_neuronal_variant_readout.tsv.gz')
target={}
for filename,ch,pos in [('integrated_15_candidate_review.csv','GRCh37_chr','GRCh38_pos'),('expanded_22_candidate_risk_gene_table.csv','GRCh38_chr','GRCh38_pos')]:
 for r in csv.DictReader((root/filename).open()):
  if not r[pos]:continue
  key=(str(r[ch]),int(r[pos]),frozenset([r['counted_allele'],r['other_allele']]))
  if key in target and target[key]!=r['variant']:raise ValueError('conflicting ID')
  target[key]=r['variant']
res=collections.defaultdict(list); n=0; keys=set()
with gzip.open(f,'rt') as z:
 for row in csv.DictReader(z,delimiter='\t'):
  n+=1
  parts=row['chr_pos_ref_alt'].split('-')
  if len(parts)!=4: continue
  c,p,ref,alt=parts
  try:key=(c.removeprefix('chr'),int(p),frozenset([ref,alt]))
  except ValueError:continue
  if key in target:
   d={'variant':target[key],'source_variant':row['chr_pos_ref_alt'],'tested_readout':row['NGN2_variant_has_readout'],'significant':row['NGN2_variant_is_significant'],'adj_p':row['NGN2_variant_adj_p-value'],'logFC':row['NGN2_variant_logFC'],'strand':row['strand'],'associated_gene':row['associated_gene_name']}
   if d not in res[target[key]]:res[target[key]].append(d)
out={'source_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'candidate_files_sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in ['integrated_15_candidate_review.csv','expanded_22_candidate_risk_gene_table.csv']},'input_candidate_contrasts':len(target),'input_rows':n,'matches':dict(res),'limitations':'Exact GRCh38 coordinate and unordered alleles only; no effect-allele orientation or exhaustive locus/proxy test; file may be annotated one row per motif; absent entries are not negative; source is neuronal reporter, not a disease outcome.'}
(root/'stage6_salomon_neuronal_candidate_screen.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'rows':n,'candidates':len(target),'matches':dict(res)},indent=2))

import csv,gzip,json,hashlib,collections
from pathlib import Path
D=Path('outputs'); q=Path('work/salomon_neuronal_variant_readout.tsv.gz'); candidates={}
for fn,c,p in [('integrated_15_candidate_review.csv','GRCh37_chr','GRCh38_pos'),('expanded_22_candidate_risk_gene_table.csv','GRCh38_chr','GRCh38_pos')]:
 for r in csv.DictReader((D/fn).open()):
  if r[p]:
   key=r['variant']; val=(r[c],int(r[p]));assert key not in candidates or candidates[key]==val;candidates[key]=val
results=collections.defaultdict(dict)
with gzip.open(q,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  fields=r['chr_pos_ref_alt'].split('-')
  if len(fields)!=4:continue
  chrom,pos,ref,alt=fields
  try: pos=int(pos)
  except ValueError:continue
  for name,(c,p) in candidates.items():
   if c==chrom.removeprefix('chr') and 0<abs(pos-p)<=10000:
    variant=r['chr_pos_ref_alt']; result=results[name].setdefault(variant,{'distance_bp':pos-p,'source_readout':r['NGN2_variant_has_readout'],'source_significant_calls':[],'source_adjusted_p':[],'associated_genes':[]})
    for source,dest in [('NGN2_variant_is_significant','source_significant_calls'),('NGN2_variant_adj_p-value','source_adjusted_p'),('associated_gene_name','associated_genes')]:
     if r[source] not in result[dest]:result[dest].append(r[source])
out={'source_sha256':hashlib.sha256(q.read_bytes()).hexdigest(),'candidate_count':len(candidates),'window_bp':10000,'matches':dict(results),'limitations':'Physical distance only, no LD or effect alignment. Source variants differ from candidates; reporter significance is not disease association, and unmeasured variants are not negative.'}
(D/'stage6_salomon_candidate_nearby_screen.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'candidates_with_nearby':len(results),'results':dict(results)},indent=2))

#!/usr/bin/env python3
"""Generate a public-coordinate request table; not evidence of cohort coverage."""
import csv,hashlib,json
from pathlib import Path
from pyliftover import LiftOver
out=Path('outputs'); forward=Path('work/hg19ToHg38.over.chain.gz'); reverse=Path('work/hg38ToHg19.over.chain.gz')
f=LiftOver(str(forward)); b=LiftOver(str(reverse))
specs=[('15_weighted',out/'unvalidated_15_variant_score_prototype.csv','chr','position_GRCh37_inferred'),('22_count',out/'unvalidated_22_shared_direction_score.csv','chr_GRCh37_assumed','position_GRCh37_assumed'),('22_weighted',out/'unvalidated_22_variant_score_prototype.csv','chr','position_GRCh37_assumed')]
rows=[]
for name,path,chrcol,poscol in specs:
 with path.open(newline='') as handle:
  for r in csv.DictReader(handle):
   chrom='chr'+r[chrcol]; p=int(r[poscol]); hits=f.convert_coordinate(chrom,p-1)
   hit=hits[0] if len(hits)==1 else None
   back=b.convert_coordinate(hit[0],hit[1]) if hit else []
   reciprocal=bool(hit and len(back)==1 and back[0][0]==chrom and back[0][1]==p-1)
   compl=str.maketrans('ACGT','TGCA'); ca=r['counted_allele']; oa=r['other_allele']
   if hit and hit[2]=='-': ca,oa=ca.translate(compl),oa.translate(compl)
   rows.append({'score':name,'variant':r['variant'],'chr37':r[chrcol],'pos37':p,'counted_allele_37':r['counted_allele'],'other_allele_37':r['other_allele'],'hg38_mapping_count':len(hits),'chr38':hit[0].removeprefix('chr') if hit else '', 'pos38':hit[1]+1 if hit else '', 'strand':hit[2] if hit else '', 'counted_allele_38_candidate':ca if hit else '', 'other_allele_38_candidate':oa if hit else '', 'unique_reciprocal_position':reciprocal,'status':'coordinate and allele candidate only; actual target variant and QC not established'})
with (out/'ncds_score_variant_request.csv').open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
summary={'rows':len(rows),'unique_rsids':len({r['variant'] for r in rows}),'single_lift':sum(r['hg38_mapping_count']==1 for r in rows),'unique_reciprocal':sum(r['unique_reciprocal_position'] for r in rows),'minus_strand':sum(r['strand']=='-' for r in rows),'file_sha256':{str(p):sha(p) for p in [forward,reverse,*[x[1] for x in specs]]},'disclaimer':'Public coordinate conversion only; not actual NCDS data, imputation QC, allele identity, or score validation.'}
(out/'ncds_score_variant_request.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

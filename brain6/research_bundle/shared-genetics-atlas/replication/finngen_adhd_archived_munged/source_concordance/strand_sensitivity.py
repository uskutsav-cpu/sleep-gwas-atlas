#!/usr/bin/env python3
"""Post-hoc strand-complement audit of the two strict allele-pair failures."""
import csv,hashlib,json,math,re
from pathlib import Path
p=Path(__file__).resolve().parent
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
sample={r['SNP']:r for r in read(p/'sample.tsv')};raw=read(p/'raw_matches.tsv');strict=read(p/'comparison.tsv'); comp={'A':'T','T':'A','C':'G','G':'C'};out=[]
for v in strict:
 if v['status']!='missing_match':continue
 x=sample[v['SNP']];a={x['A1'],x['A2']};assert a not in ({'A','T'},{'C','G'}) and a.issubset(comp)
 matches=[r for r in raw if v['SNP'] in re.split(r'[,;]',r['rsids']) and {comp[x['A1']],comp[x['A2']]}=={r['ref'],r['alt']}]
 row={'SNP':v['SNP'],'processed_A1':x['A1'],'processed_A2':x['A2'],'processed_Z':x['Z'],'status':'unresolved','raw_chr':'','raw_pos38':'','raw_ref':'','raw_alt':'','aligned_beta_over_se':'','difference':''}
 if len(matches)==1:
  r=matches[0];z=float(r['beta'])/float(r['sebeta'])*(1 if comp[x['A1']]==r['alt'] else -1)
  row.update(status='complement_match',raw_chr=r['#chrom'],raw_pos38=r['pos'],raw_ref=r['ref'],raw_alt=r['alt'],aligned_beta_over_se=z,difference=float(x['Z'])-z)
 out.append(row)
with (p/'strand_sensitivity.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
receipt={'selected_from':'strict missing_match rows after strict comparison','number':len(out),'complement_matches':sum(r['status']=='complement_match' for r in out),'maximum_absolute_z_difference':max(abs(float(r['difference'])) for r in out if r['status']=='complement_match'),'input_sha256':{f:hashlib.sha256((p/f).read_bytes()).hexdigest() for f in ('sample.tsv','raw_matches.tsv','comparison.tsv')},'limitation':'Coordinate-free rsID and strand compatibility do not prove molecular identity or original processing lineage.'}
(p/'strand_sensitivity_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

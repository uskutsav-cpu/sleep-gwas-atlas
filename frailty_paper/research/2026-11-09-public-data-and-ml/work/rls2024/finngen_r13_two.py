import json,hashlib,csv
from pathlib import Path
import pysam
u='https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_G6_RLS.gz';t=pysam.TabixFile(u);header=list(t.header);assert len(header)==1;fields=header[0].lstrip('#').split('\t');assert fields==['chrom','pos','ref','alt','rsids','nearest_genes','pval','mlogp','beta','sebeta','af_alt','af_alt_cases','af_alt_controls']
out=[]
for rs,pos,ref,alt in [('rs113851554',66523432,'G','T'),('rs78206187',66795102,'A','G')]:
 hits=[]
 for line in t.fetch('2',pos-1,pos):
  vals=line.split('\t');assert len(vals)==len(fields);r=dict(zip(fields,vals))
  if int(r['pos'])==pos and r['ref']==ref and r['alt']==alt:hits.append(r)
 assert len(hits)==1 and rs in hits[0]['rsids'].split(',');r=hits[0];out.append({'candidate':rs,'chr':r['chrom'],'position':r['pos'],'ref':ref,'alt':alt,'beta_alt':r['beta'],'se':r['sebeta'],'p':r['pval'],'af_alt':r['af_alt']})
p=Path('outputs/stage6_finngen_r13_rls_two.csv')
with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
r={'url':u,'header':header[0],'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'method':'Remote tabix, exact GRCh38 ref/alt and rsID','full_file_checksum_verified':False,'limitations':'Postselected follow-up, phenotype differs and participant overlap unknown; target cohort is not an independent score-validation dataset.'}
Path('outputs/stage6_finngen_r13_rls_two_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(out,indent=2))

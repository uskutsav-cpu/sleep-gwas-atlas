import csv,json,hashlib
from pathlib import Path
import pysam
sources={'Akcimen':'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90435001-GCST90436000/GCST90435387/harmonised/GCST90435387.h.tsv.gz','FinnGen':'https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_G6_RLS.gz'}
targets=[('rs187012067',66511859,'T','C'),('rs62144055',66522178,'A','C'),('rs113851554',66523432,'G','T'),('rs12469063',66537176,'A','G'),('rs78206187',66795102,'A','G')]
fields={'Akcimen':['chrom','pos','ea','oa','beta','se','eaf','p','rsid','unused','hm_coordinate_conversion','hm_code','variant_id'],'FinnGen':['chrom','pos','ref','alt','rsid','nearest_genes','p','mlogp','beta','se','af_alt','af_alt_cases','af_alt_controls']};out=[]
for name,u in sources.items():
 t=pysam.TabixFile(u)
 for rs,pos,ref,alt in targets:
  hit=[]
  for line in t.fetch('2',pos-1,pos):
   f=line.rstrip('\r\n').split('\t');assert len(f)==len(fields[name]);r=dict(zip(fields[name],f))
   if name=='Akcimen':
    if r['variant_id']!=f'2_{pos}_{ref}_{alt}' or {r['ea'],r['oa']}!={ref,alt}:continue
    af=float(r['eaf']) if r['ea']==alt else 1-float(r['eaf']);beta=float(r['beta'])*(1 if r['ea']==alt else -1);code=r['hm_code']
   else:
    if r['ref']!=ref or r['alt']!=alt:continue
    af=float(r['af_alt']);beta=float(r['beta']);code='not_applicable'
   if rs not in r['rsid'].split(','):continue
   assert 0<=af<=1;hit.append({'study':name,'variant':rs,'position_GRCh38':pos,'ref':ref,'alt':alt,'beta_alt':beta,'af_alt':af,'harmonization_code':code,'p':float(r['p'])})
  assert len(hit)==1,(name,rs,len(hit));out+=hit
p=Path('outputs/stage6_rls_region_allele_frequency_spotcheck.csv')
with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
Path('outputs/stage6_rls_region_allele_frequency_spotcheck_receipt.json').write_text(json.dumps({'sources':sources,'selected_variants':len(targets),'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'full_files_checksum_verified':False,'caveat':'Cross-population allele-frequency agreement is a plausibility check only; cannot prove coding, allele-specific causality, or explanation of association discordance.'},indent=2)+'\n')
print(json.dumps(out,indent=2))

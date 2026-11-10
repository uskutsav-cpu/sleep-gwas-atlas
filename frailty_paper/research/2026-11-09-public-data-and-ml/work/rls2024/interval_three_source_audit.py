import csv,json,hashlib,math,collections
from pathlib import Path
import pysam
src=Path('outputs/stage6_rls_region_source_concordance.csv');target=[r for r in csv.DictReader(src.open()) if r['both_gws']=='1'];assert len(target)==42
u='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90399001-GCST90400000/GCST90399571/harmonised/GCST90399571.h.tsv.gz';hdr=['chrom','pos','ea','oa','beta','se','eaf','p','rsid','hm_coordinate_conversion','hm_code','variant_id'];t=pysam.TabixFile(u);data=collections.defaultdict(list)
for line in t.fetch('2',66499999,66820000):
 f=line.rstrip('\r\n').split('\t');assert len(f)==len(hdr);r=dict(zip(hdr,f))
 if r['hm_code'] not in ('10','11','12','13','5','6'):continue
 try:chrom,pos,ref,alt=r['variant_id'].split('_');pos=int(pos);p=float(r['p']);beta=float(r['beta']);se=float(r['se'])
 except (ValueError,TypeError):continue
 if chrom!='2' or pos!=int(r['pos']) or {ref,alt}!={r['ea'],r['oa']} or not math.isfinite(beta) or not 0<=p<=1 or not math.isfinite(se) or se<=0:continue
 data[(pos,ref,alt)].append((beta*(1 if r['ea']==alt else -1),se,p,r['rsid'],r['hm_code']))
out=[]
for row in target:
 k=(int(row['pos']),row['ref'],row['alt']);hits=data.get(k,[]);o={'pos':row['pos'],'ref':row['ref'],'alt':row['alt'],'akcimen_rsid':row['akcimen_rsid'],'akcimen_beta_alt':row['akcimen_beta_alt'],'finngen_beta_alt':row['finngen_beta_alt'],'first_two_same_sign':row['same_sign'],'interval_status':'missing_or_duplicate'}
 if len(hits)==1:
  beta,se,p,rs,hm=hits[0];o.update(interval_status='matched',interval_rsid=rs,interval_beta_alt=beta,interval_se=se,interval_p=p,interval_hm_code=hm,interval_agrees_akcimen=int(beta*float(row['akcimen_beta_alt'])>0),interval_agrees_finngen=int(beta*float(row['finngen_beta_alt'])>0))
 out.append(o)
fields=['pos','ref','alt','akcimen_rsid','akcimen_beta_alt','finngen_beta_alt','first_two_same_sign','interval_status','interval_rsid','interval_beta_alt','interval_se','interval_p','interval_hm_code','interval_agrees_akcimen','interval_agrees_finngen'];p=Path('outputs/stage6_rls_three_study_region.csv')
with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
sub=[r for r in out if r['first_two_same_sign']=='0'];m=[r for r in sub if r['interval_status']=='matched'];nom=[r for r in m if r['interval_p']<.05]
receipt={'source':u,'selected_rows':len(out),'first_two_discordant':len(sub),'discordant_interval_matched':len(m),'discordant_interval_nominal_p_lt_05':len(nom),'nominal_interval_agrees_akcimen':sum(r['interval_agrees_akcimen'] for r in nom),'nominal_interval_agrees_finngen':sum(r['interval_agrees_finngen'] for r in nom),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'full_file_checksum_verified':False,'limitation':'Variants selected by significance in first two studies and in a postselected locus, with correlated rows; study overlap unknown and nominal p values do not provide independent hypothesis tests.'};Path('outputs/stage6_rls_three_study_region_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

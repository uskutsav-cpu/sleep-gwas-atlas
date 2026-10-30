#!/usr/bin/env python3
"""Stream-check prototype variants against user project's existing harmonized GWAS files.
No individual-level genotype data are read. Paths point to the user's existing external drive.
"""
import argparse,csv,gzip,hashlib,json
from pathlib import Path
import pandas as pd
p=argparse.ArgumentParser();p.add_argument('--workspace',type=Path,default=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1'));args=p.parse_args();P=Path(__file__).resolve().parent
s=pd.concat([pd.read_csv(P/f'score_{n}_variants_standardized.csv').assign(score=n) for n in (15,22)],ignore_index=True);ids=set(s.rsid);data={};receipts=[]
paths={'FI':args.workspace/'analysis-workspace/frailty_v1/harmonized/frailty.harmonized.tsv.gz','insomnia':args.workspace/'analysis-workspace/frailty_v1/harmonized/insomnia.harmonized.tsv.gz','FFS':args.workspace/'frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv'}
for trait,path in paths.items():
 if not path.exists():receipts.append({'trait':trait,'path':str(path),'status':'missing'});continue
 h=hashlib.sha256();
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(4*1024*1024),b''):h.update(chunk)
 op=gzip.open if path.suffix=='.gz' else open;hits=[];scanned=0
 with op(path,'rt',encoding='utf8',newline='') as f:
  rd=csv.DictReader(f,delimiter='\t');needed={'SNP','A1','A2','BETA'};assert needed<=set(rd.fieldnames),(trait,rd.fieldnames)
  for row in rd:
   scanned+=1
   if row['SNP'] in ids:hits.append(row)
 data[trait]=hits;receipts.append({'trait':trait,'path':str(path),'sha256':h.hexdigest(),'n_rows':scanned,'n_matching_rsid_rows':len(hits),'status':'read','interpretation':'project GWAS summary association data; A1 effect convention is source-dependent'})
rows=[]
for _,r in s.iterrows():
 for trait,field in [('FI','FI_beta_weight'),('insomnia','insomnia_beta_weight'),('FFS','FFS_beta_weight_inferred')]:
  source=data.get(trait);hits=[] if source is None else [x for x in source if x['SNP']==r.rsid and str(x.get('CHR',x.get('CHROM')))==str(r.chromosome) and str(x.get('BP',x.get('POS')))==str(r.position)]
  o={'score':int(r.score),'rsid':r.rsid,'trait':trait,'model_counted_allele':r.effect_counted_allele,'model_other_allele':r.other_allele,'model_weight':r[field],'source_file':str(paths[trait]),'exact_rsid_position_rows':len(hits),'status':'source_missing' if source is None else 'not_found' if not hits else 'multiple' if len(hits)>1 else 'checked'}
  if len(hits)==1:
   h=hits[0];a,b=h['A1'].upper(),h['A2'].upper();align='same' if (a,b)==(r.effect_counted_allele,r.other_allele) else 'swapped' if (a,b)==(r.other_allele,r.effect_counted_allele) else 'mismatch';beta=float(h['BETA']);aligned=beta if align=='same' else -beta if align=='swapped' else None
   o.update(source_A1=a,source_A2=b,source_beta=beta,source_p=h.get('P',''),alignment=align,aligned_beta=aligned,absolute_difference=abs(aligned-r[field]) if aligned is not None else None,weight_agrees=abs(aligned-r[field])<1e-10 if aligned is not None else False)
  rows.append(o)
out=pd.DataFrame(rows);out.to_csv(P/'supplement_gwas_weight_audit.csv',index=False);(P/'gwas_summary_source_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
print(out.groupby(['trait','status']).size().to_string());print('weight agreement:',out.groupby('trait')['weight_agrees'].sum().to_dict());print('Note: GWAS file QC filters and A1 semantics limit interpretation')

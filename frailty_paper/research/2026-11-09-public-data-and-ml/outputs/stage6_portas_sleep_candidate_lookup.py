#!/usr/bin/env python3
"""Remote indexed exact-allele lookups; no whole-file checksum verification."""
import csv,json,hashlib,concurrent.futures,sys
from pathlib import Path
import pandas as pd,pysam
O=Path('outputs');W=Path('work'); d=pd.read_csv(O/'candidate_finngen_hfrs_full_gwas_lookup.csv'); inputs=[{'variant':r.variant,'chr':str(r.GRCh38_chr),'pos':int(r.GRCh38_pos),'a':r.candidate_counted_allele,'b':r.candidate_other_allele} for r in d.itertuples()];inputs.extend([{'variant':'rs35624276','chr':'7','pos':114400048,'a':'G','b':'A'},{'variant':'rs1799252','chr':'7','pos':114227537,'a':'A','b':'T'}]);
# The signal lead position is independently checked from gnomAD's GRCh37 coordinate converted to GRCh38.
from pyliftover import LiftOver
l=LiftOver(str(W/'hg19ToHg38.over.chain.gz'));m=l.convert_coordinate('chr7',114040102); assert len(m)==1 and m[0][0]=='chr7' and m[0][2]=='+';inputs[-2]['pos']=m[0][1]+1
accs=['GCST90824106','GCST90824107'];base='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90824001-GCST90825000/'
def job(acc):
 u=base+acc+'/harmonised/'+acc+'.h.tsv.gz';p=W/(acc+'.h.tsv.gz.tbi');f=pysam.TabixFile(u,index=str(p));out=[]
 for q in inputs:
  lines=list(f.fetch(q['chr'],q['pos']-1,q['pos'])); matches=[]
  for line in lines:
   c=line.split('\t');assert len(c)>=8
   if c[1]==str(q['pos']) and {c[2],c[3]}=={q['a'],q['b']}:matches.append(c)
  x={'accession':acc,'variant':q['variant'],'chr38':q['chr'],'pos38':q['pos'],'counted_allele':q['a'],'other_allele':q['b'],'matched_allele_rows':len(matches),'position_rows':len(lines)}
  if len(matches)==1:
   c=matches[0];x.update(source_effect_allele=c[2],source_other_allele=c[3],beta_counted=float(c[4])*(1 if c[2]==q['a'] else -1),standard_error=float(c[5]),source_eaf=float(c[6]),p_value=float(c[7]),extra_columns='|'.join(c[8:]))
  out.append(x);print(acc,q['variant'],len(matches),flush=True)
 return out
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as e:result=sum(list(e.map(job,accs)),[])
pd.DataFrame(result).to_csv(O/'stage6_portas_sleep_candidate_lookup.csv',index=False);(O/'stage6_portas_sleep_candidate_receipt.json').write_text(json.dumps({'traits':{'GCST90824106':'REM duration','GCST90824107':'NREM duration'},'source':base,'reference_genome':'GRCh38 harmonised','indices_md5':{a:hashlib.md5((W/(a+'.h.tsv.gz.tbi')).read_bytes()).hexdigest() for a in accs},'target_variants':len(inputs),'whole_file_checksum_verified':False,'matched_by':'GRCh38 position and exact ordered allele pair (order-independent), not rsID','limitation':'sex-combined UKB; overlaps possible; missing variants do not mean null; SNP identity at multiallelic loci not fully resolved'},indent=2)+'\n');print('done',len(result),flush=True)

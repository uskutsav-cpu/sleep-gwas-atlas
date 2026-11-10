from docx import Document
from pathlib import Path
import csv,hashlib,collections
p=Path('work/song_supplement.docx'); table=Document(p).tables[3]
rows=[]
for i,row in enumerate(table.rows[1:],start=2):
 x=[c.text.strip() for c in row.cells]; assert len(x)==8 and x[0].startswith('rs')
 snp,chrom,pos,ea,oa,pfi,pins,pc=x
 rows.append({'source_row':i,'variant_id':snp,'chromosome_reported':int(chrom),'position_reported':int(pos),'genome_build':'not resolved from inspected table','effect_allele_reported':ea,'other_allele_reported':oa,'mtag_fi_p':pfi,'mtag_insomnia_p':pins,'cpassoc_p':pc,'approx_chr6_MHC_flag':'yes' if int(chrom)==6 and 25000000<=int(pos)<=34000000 else 'no','effect_direction':'not in source table','FFS_test':'not in source table','valid_score_weight':'no, no effect weights','interpretation':'Published exploratory cross-trait candidate; correlated MTAG statistics and CPASSOC do not prove same causal variant or direction; LD independence unverified','source':'https://doi.org/10.3389/fnagi.2024.1358996','source_supplement':'Data_Sheet_1.docx Table S4'})
assert len(rows)==124 and len({x['variant_id'] for x in rows})==124
assert all(float(x[k])<5e-8 for x in rows for k in ('mtag_fi_p','mtag_insomnia_p','cpassoc_p'))
rows.sort(key=lambda x:(float(x['cpassoc_p']),x['variant_id']))
for i,x in enumerate(rows,1):x['rank_by_cpassoc_p']=i
for name,sub in [('outputs/published_insomnia_fi_124_cross_trait_candidates.csv',rows),('outputs/published_insomnia_fi_top30_candidates.csv',rows[:30])]:
 with open(name,'w',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=rows[0].keys());wr.writeheader();wr.writerows(sub)
print('source docx sha256',hashlib.sha256(p.read_bytes()).hexdigest())
print('candidates',len(rows),'top30 MHC_flag',sum(x['approx_chr6_MHC_flag']=='yes' for x in rows[:30]))
print('top30 distinct chromosomes',len(set(x['chromosome_reported'] for x in rows[:30])))
print('known row ranks',[(x['variant_id'],x['rank_by_cpassoc_p']) for x in rows if x['variant_id'] in {'rs34290943','rs10865954'}])

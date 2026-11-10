import csv,gzip,hashlib,json,pathlib
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scipy.stats import t
p=pathlib.Path('work/geo_insomnia/GSE208668_series_matrix.txt.gz')
target=['FOXP2','SMAD5','RBM6','MEIS1']; meta={};values={};active=False
with gzip.open(p,'rt',newline='') as f:
 for row in csv.reader(f,delimiter='\t'):
  if not row: continue
  if row[0]=='!series_matrix_table_begin': active=True;continue
  if row[0]=='!series_matrix_table_end':break
  if row[0]=='!Sample_characteristics_ch1':
   field=row[1].split(': ',1)[0];meta[field]=[v.split(': ',1)[1] for v in row[1:]]
  if active and row[0] in target:
   if row[0] in values:raise ValueError('duplicate gene symbol')
   values[row[0]]=np.array(row[1:],float)
assert set(values)==set(target) and len(meta['age'])==42
ins=np.array([x=='yes' for x in meta['insomnia']],float)
age=np.array(meta['age'],float); sex=np.array([x=='male' for x in meta['gender']],float)
assert sum(ins)==17 and sum(1-ins)==25
x=np.column_stack([np.ones(len(ins)),ins,age,sex]); rows=[]
for gene in target:
 fit=sm.OLS(values[gene],x).fit(cov_type='HC3',use_t=True)
 rows.append({'gene':gene,'n':len(ins),'insomnia_minus_control_log2':fit.params[1],'hc3_se':fit.bse[1],'p':fit.pvalues[1],'ci_low':fit.conf_int()[1,0],'ci_high':fit.conf_int()[1,1]})
for row,q in zip(rows,multipletests([r['p'] for r in rows],method='fdr_bh')[1]):row['bh_q_four']=q
with open('outputs/stage6_gse208668_candidate_expression.csv','w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
receipt={'source':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208668/matrix/GSE208668_series_matrix.txt.gz','source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'n':len(ins),'cases':int(sum(ins)),'controls':int(sum(1-ins)),'genes':target,'model':'OLS, insomnia + age + male; HC3 with residual t degrees of freedom','data':'submitted log2 gene-symbol matrix','warning':'Insomnia group and source study perfectly confounded; not causality, genotype evidence or score validation.'}
pathlib.Path('outputs/stage6_gse208668_candidate_expression.json').write_text(json.dumps(receipt,indent=2)+'\n')
for r in rows: print(r)

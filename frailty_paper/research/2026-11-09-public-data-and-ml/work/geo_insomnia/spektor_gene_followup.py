#!/usr/bin/env python3
"""Post hoc whole-blood-selected panel in public PBMC data; not causal replication."""
import csv,gzip,hashlib,json,pathlib
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
with open('work/spektor2026/significant_genes.csv',newline='') as f:source={r['gene']:r for r in csv.DictReader(f)}
assert len(source)==244
p=pathlib.Path('work/geo_insomnia/GSE208668_series_matrix.txt.gz'); meta={}; values={};active=False
with gzip.open(p,'rt',newline='') as f:
 for row in csv.reader(f,delimiter='\t'):
  if not row:continue
  if row[0]=='!series_matrix_table_begin':active=True;continue
  if row[0]=='!series_matrix_table_end':break
  if row[0]=='!Sample_characteristics_ch1':
   field=row[1].split(': ',1)[0];meta[field]=[v.split(': ',1)[1] for v in row[1:]]
  if active and row[0] in source:
   if row[0] in values:raise ValueError('duplicate '+row[0])
   vals=np.array(row[1:],float)
   if not np.isfinite(vals).all():raise ValueError('nonfinite '+row[0])
   values[row[0]]=vals
ins=np.array([x=='yes' for x in meta['insomnia']],float);age=np.array(meta['age'],float);sex=np.array([x=='male' for x in meta['gender']],float)
assert len(ins)==42 and sum(ins)==17
x=np.column_stack([np.ones(len(ins)),ins,age,sex]);rows=[]
for gene in sorted(values):
 fit=sm.OLS(values[gene],x).fit(cov_type='HC3',use_t=True); v=source[gene]
 rows.append({'gene':gene,'source_log2fc':float(v['source_log2fc']),'source_padj':float(v['source_padj']),'pbmc_log2_difference':fit.params[1],'pbmc_hc3_se':fit.bse[1],'pbmc_p':fit.pvalues[1],'pbmc_ci_low':fit.conf_int()[1,0],'pbmc_ci_high':fit.conf_int()[1,1],'direction_same':np.sign(fit.params[1])==np.sign(float(v['source_log2fc']))})
for r,q in zip(rows,multipletests([r['pbmc_p'] for r in rows],method='fdr_bh')[1]):r['pbmc_bh_q_in_matched_panel']=q
with open('outputs/stage6_spektor_gse208668_gene_comparison.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
selected=[r for r in rows if r['pbmc_bh_q_in_matched_panel']<.05]
receipt={'source_panel':244,'genes_matched':len(rows),'genes_unmatched':len(source)-len(rows),'significant_pbmc_bh':len(selected),'significant_pbmc_same_direction':int(sum(r['direction_same'] for r in selected)),'all_matched_same_direction':int(sum(r['direction_same'] for r in rows)),'matrix_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'design':'Submitted log2 PBMC expression, group+age+sex OLS HC3; recruitment cohort perfectly confounded with group','limitations':'Post hoc selection, whole blood vs PBMC, independent participant status unknown, no genotype or frailty outcome.'}
pathlib.Path('outputs/stage6_spektor_gse208668_gene_comparison.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

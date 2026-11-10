#!/usr/bin/env python3
"""Post hoc screen of published insomnia-expression genes in paired wakefulness data."""
import csv,gzip,hashlib,json,re
import numpy as np
from scipy.stats import ttest_1samp,t
from statsmodels.stats.multitest import multipletests
with open('work/spektor2026/significant_genes.csv',newline='') as f: source={r['gene']:r for r in csv.DictReader(f)}
assert len(source)==244
P='work/geo_sleep_wake/GSE37667_series_matrix.txt.gz';A='work/geo_sleep_wake/GPL570.annot.gz';probe_to_genes={};active=False
with gzip.open(A,'rt',errors='replace') as f:
 for row in csv.reader(f,delimiter='\t'):
  if row and row[0]=='!platform_table_begin':active=True;continue
  if row and row[0]=='!platform_table_end':break
  if active and len(row)>2 and row[2] in source:probe_to_genes.setdefault(row[0],set()).add(row[2])
values={g:[] for g in source};titles=None;active=False
with gzip.open(P,'rt',newline='') as f:
 for row in csv.reader(f,delimiter='\t'):
  if not row:continue
  if row[0]=='!Sample_title':titles=row[1:]
  if row[0]=='!series_matrix_table_begin':active=True;continue
  if row[0]=='!series_matrix_table_end':break
  if active and row[0] in probe_to_genes:
   y=np.asarray(row[1:],float)
   if len(y)!=27 or not np.isfinite(y).all():raise ValueError('invalid values')
   for gene in probe_to_genes[row[0]]:values[gene].append((row[0],y))
assert len(titles)==27
ix={}
for i,title in enumerate(titles):
 m=re.fullmatch(r'(baseline|Prolonged_wakefulness|Sleep recovery)_volunteer([1-9])',title)
 if not m:raise ValueError(title)
 key=(m.group(1),m.group(2));assert key not in ix;ix[key]=i
assert len(ix)==27
rows=[]
for g in sorted(source):
 if not values[g]:continue
 probe,y=max(values[g],key=lambda z:np.mean(z[1]));d=np.array([y[ix[('Prolonged_wakefulness',str(k))]]-y[ix[('baseline',str(k))]] for k in range(1,10)]);se=d.std(ddof=1)/3;mean=d.mean();p=float(ttest_1samp(d,0).pvalue)
 if not np.isfinite(p): raise ValueError('nonfinite p '+g)
 rows.append({'gene':g,'probe':probe,'eligible_probes':len(values[g]),'n_pairs':9,'source_insomnia_log2fc':float(source[g]['source_log2fc']),'wake_minus_baseline_log2':float(mean),'se':float(se),'ci_low':float(mean-t.ppf(.975,8)*se),'ci_high':float(mean+t.ppf(.975,8)*se),'p':p,'same_direction':bool(np.sign(mean)==np.sign(float(source[g]['source_log2fc'])))})
assert rows
for r,q in zip(rows,multipletests([r['p'] for r in rows],method='fdr_bh')[1]):r['bh_q']=float(q)
with open('outputs/stage6_spektor_gse37667_gene_comparison.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
hits=[r for r in rows if r['bh_q']<.05]
receipt={'source_genes':len(source),'genes_with_probe':len(rows),'unmatched':len(source)-len(rows),'significant_bh':len(hits),'significant_same_direction':sum(r['same_direction'] for r in hits),'all_same_direction':sum(r['same_direction'] for r in rows),'n_pairs':9,'probe_rule':'single exact gene symbol and maximal mean over all 27 samples','series_sha256':hashlib.sha256(open(P,'rb').read()).hexdigest(),'annotation_sha256':hashlib.sha256(open(A,'rb').read()).hexdigest()}
open('outputs/stage6_spektor_gse37667_gene_comparison.json','w').write(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

import csv,gzip,hashlib,json,re
import numpy as np
from scipy.stats import ttest_1samp,t
from statsmodels.stats.multitest import multipletests
P='work/geo_sleep_wake/GSE37667_series_matrix.txt.gz';A='work/geo_sleep_wake/GPL570.annot.gz'; genes=['FOXP2','SMAD5','RBM6','MEIS1']; probes={g:set() for g in genes}
with gzip.open(A,'rt',errors='replace') as f:
 active=False
 for row in csv.reader(f,delimiter='\t'):
  if row and row[0]=='!platform_table_begin':active=True;continue
  if row and row[0]=='!platform_table_end':break
  if active and len(row)>2 and row[2] in probes:probes[row[2]].add(row[0])
values={g:[] for g in genes};titles=None;active=False
with gzip.open(P,'rt',newline='') as f:
 for row in csv.reader(f,delimiter='\t'):
  if not row:continue
  if row[0]=='!Sample_title':titles=row[1:]
  if row[0]=='!series_matrix_table_begin':active=True;continue
  if row[0]=='!series_matrix_table_end':break
  if active:
   for g in genes:
    if row[0] in probes[g]:values[g].append((row[0],np.asarray(row[1:],float)))
assert len(titles)==27
ix={}
for i,title in enumerate(titles):
 m=re.fullmatch(r'(baseline|Prolonged_wakefulness|Sleep recovery)_volunteer([1-9])',title)
 if not m:raise ValueError(title)
 key=(m.group(1),m.group(2));assert key not in ix;ix[key]=i
assert len(ix)==27
rows=[]
for g in genes:
 assert values[g]
 name,y=max(values[g],key=lambda z:np.mean(z[1]))
 diffs=np.array([y[ix[('Prolonged_wakefulness',str(k))]]-y[ix[('baseline',str(k))]] for k in range(1,10)])
 rec=np.array([y[ix[('Sleep recovery',str(k))]]-y[ix[('baseline',str(k))]] for k in range(1,10)])
 mean=diffs.mean();se=diffs.std(ddof=1)/3;crit=t.ppf(.975,8)
 rows.append({'gene':g,'probe':name,'n_probes_eligible':len(values[g]),'n_pairs':9,'wake_minus_baseline_log2':float(mean),'paired_se':float(se),'ci_low':float(mean-crit*se),'ci_high':float(mean+crit*se),'p':float(ttest_1samp(diffs,0).pvalue),'recovery_minus_baseline_mean_log2':float(rec.mean())})
for row,q in zip(rows,multipletests([r['p'] for r in rows],method='fdr_bh')[1]):row['bh_q_four']=float(q)
with open('outputs/stage6_gse37667_candidate_paired.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
receipt={'series_sha256':hashlib.sha256(open(P,'rb').read()).hexdigest(),'annotation_sha256':hashlib.sha256(open(A,'rb').read()).hexdigest(),'rule':'unique exact symbol, highest mean probe across all 27 samples','n_subjects':9,'analysis':'paired t test prolonged wakefulness minus baseline, recovery descriptive','warning':'postselected small sample, no randomized time control; sleep loss not clinical insomnia, no genetics or frailty.'}
open('outputs/stage6_gse37667_candidate_paired.json','w').write(json.dumps(receipt,indent=2)+'\n')
for r in rows:print(r)

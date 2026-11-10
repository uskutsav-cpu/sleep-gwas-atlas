import csv,gzip,json,math,collections
from pathlib import Path
import numpy as np
from scipy import stats
probes=json.loads(Path('outputs/stage6_gse40366_candidate_probe_map.json').read_text())['rows'];by=collections.defaultdict(list)
for r in probes:by[r['gene']].append(r['probe'])
char=[];values={}
with gzip.open('work/geo_vitality/GSE40366.txt.gz','rt') as f:
 reader=csv.reader(f,delimiter='\t')
 for r in reader:
  if not r:continue
  if r[0]=='!Sample_geo_accession':samples=r[1:]
  if r[0]=='!Sample_characteristics_ch1':char.append(r[1:])
  if r[0]=='!series_matrix_table_begin':break
 assert next(reader)[1:]==samples
 selected={r['probe'] for r in probes}
 for r in reader:
  if r and r[0]=='!series_matrix_table_end':break
  if r and r[0] in selected:assert r[0] not in values;values[r[0]]=np.array([float(x) for x in r[1:]])
meta=[{} for _ in samples]
for row in char:
 assert len(row)==len(samples)
 for i,v in enumerate(row):
  if ': ' in v:
   k,x=v.split(': ',1);meta[i][k.lower()]=x
allfrailty=[i for i,x in enumerate(meta) if x.get('frailty index') in ['non-frail','pre-frail','frail']];assert len(allfrailty)==151
idx=[i for i in allfrailty if meta[i]['frailty index'] in ['non-frail','frail']];assert len(idx)==67
frail=np.array([int(meta[i]['frailty index']=='frail') for i in idx]);assert frail.sum()==29
sexnames=[meta[i]['sex'].lower() for i in idx];assert set(sexnames)<=set(['male','female']);sex=np.array([int(x=='male') for x in sexnames]);X=np.column_stack([np.ones(len(idx)),frail,sex]);inv=np.linalg.inv(X.T@X);hat=np.sum((X@inv)*X,axis=1);out=[]
for gene in ['FOXP2','SMAD5','RBM6','MEIS1']:
 options=[p for p in by[gene] if p in values];assert options
 probe=max(options,key=lambda p:values[p][allfrailty].mean());y=values[probe][idx];beta=inv@X.T@y;err=y-X@beta;cov=inv@(X.T@((err/(1-hat))[:,None]**2*X))@inv;se=math.sqrt(cov[1,1]);pval=float(2*stats.t.sf(abs(beta[1]/se),df=len(idx)-3));out.append({'gene':gene,'probe_selected':probe,'available_probes':len(options),'adjusted_difference':float(beta[1]),'hc3_se':se,'p_unadjusted':pval})
order=sorted(range(4),key=lambda i:out[i]['p_unadjusted']);p=1
for k in range(3,-1,-1):
 i=order[k];p=min(p,out[i]['p_unadjusted']*4/(k+1));out[i]['bh_q_four_genes']=p
Path('outputs/stage6_gse40366_frailty_expression.json').write_text(json.dumps({'results':out,'n_frail':29,'n_nonfrail':38,'n_prefrail_excluded':84,'probe_selection':'highest mean expression across all 151 nonagenarians with available categories, ignoring group','method':'GEO quantile-normalized log2 expression; OLS group and sex, HC3 standard errors, BH over four preselected project genes','limits':'Age all nonagenarian without exact age adjustment; survivor selection, blood cell composition, batch and smoking unadjusted, unknown age covariates and case-selection; no genotype linkage or insomnia outcomes; not score validation'},indent=2)+'\n');print(json.dumps(out,indent=2))

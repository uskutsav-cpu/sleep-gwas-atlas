import csv, gzip, io, json, math, collections
from pathlib import Path
import numpy as np
import subprocess
from scipy import stats
from statsmodels.stats.multitest import multipletests
source=list(csv.DictReader(open('work/spektor2026/significant_genes.csv')))
assert len(source)==244 and len({r['gene'] for r in source})==244
targets={r['gene'] for r in source}; by=collections.defaultdict(list)
# Only a historical approved symbol which uniquely maps to one approved HGNC gene
# can be used to rescue a missing exact symbol; ambiguous and current-symbol collisions excluded.
hgnc=list(csv.DictReader(open('work/geo_vitality/hgnc_complete_set.txt'),delimiter='\t'))
approved={r['symbol']:r for r in hgnc if r['status']=='Approved'}
inverse=collections.defaultdict(set)
for g,r in approved.items():
 for v in filter(None, (r['prev_symbol']+'|'+r['alias_symbol']).split('|')):inverse[v].add(g)
for g in approved:inverse[g].add(g)
alias_map={}
for g in targets:
 if g in approved:
  for prev in filter(None,approved[g]['prev_symbol'].split('|')):
   if inverse[prev]=={g} and prev not in targets:alias_map[prev]=g

u='https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL15nnn/GPL15988/soft/GPL15988_family.soft.gz'
proc=subprocess.Popen(['curl','-fsSL',u],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
with gzip.GzipFile(fileobj=proc.stdout) as z:
 f=io.TextIOWrapper(z)
 for line in f:
  if line.startswith('!platform_table_begin'):break
 header=next(f).rstrip('\n').split('\t');assert header[:4]==['ID','GB_ACC','GENE SYMBOL','GENE DESCRIPTION']
 for line in f:
  if line.startswith('!platform_table_end'):break
  x=line.rstrip('\n').split('\t')
  if len(x)>=3 and x[2] in targets:by[x[2]].append((x[0],x[2],'exact'))
  elif len(x)>=3 and x[2] in alias_map:by[alias_map[x[2]]].append((x[0],x[2],'previous_symbol'))
proc.terminate();proc.wait()
char=[];values={};selected={p for ps in by.values() for p,_,_ in ps}
with gzip.open('work/geo_vitality/GSE40366.txt.gz','rt') as f:
 reader=csv.reader(f,delimiter='\t')
 for row in reader:
  if not row:continue
  if row[0]=='!Sample_geo_accession':samples=row[1:]
  if row[0]=='!Sample_characteristics_ch1':char.append(row[1:])
  if row[0]=='!series_matrix_table_begin':break
 assert next(reader)[1:]==samples
 for row in reader:
  if row and row[0]=='!series_matrix_table_end':break
  if row and row[0] in selected:
   assert row[0] not in values
   values[row[0]]=np.array([float(x) for x in row[1:]])
meta=[{} for _ in samples]
for row in char:
 assert len(row)==len(samples)
 for i,v in enumerate(row):
  if ': ' in v:
   k,x=v.split(': ',1);meta[i][k.lower()]=x
allidx=[i for i,x in enumerate(meta) if x.get('frailty index') in ['non-frail','pre-frail','frail']];assert len(allidx)==151
idx=[i for i in allidx if meta[i]['frailty index'] in ['non-frail','frail']];assert len(idx)==67
frail=np.array([int(meta[i]['frailty index']=='frail') for i in idx]);assert frail.sum()==29
sex=np.array([int(meta[i]['sex'].lower()=='male') for i in idx]);assert all(meta[i]['sex'].lower() in ('male','female') for i in idx)
X=np.column_stack([np.ones(len(idx)),frail,sex]);inv=np.linalg.inv(X.T@X);hat=np.sum((X@inv)*X,axis=1);out=[]
for row in source:
 options=[(p,sym,kind) for p,sym,kind in by[row['gene']] if p in values]
 if not options:continue
 # Preserve exact-symbol analyses when exact probes exist; use previous symbols only as a rescue.
 exact=[v for v in options if v[2]=='exact']
 if exact: options=exact
 probe,platform_symbol,mapping=max(options,key=lambda v:values[v[0]][allidx].mean());y=values[probe][idx];beta=inv@X.T@y;err=y-X@beta;cov=inv@(X.T@((err/(1-hat))[:,None]**2*X))@inv;se=math.sqrt(cov[1,1]);p=float(2*stats.t.sf(abs(beta[1]/se),df=len(idx)-3))
 out.append(dict(gene=row['gene'],source_insomnia_log2fc=float(row['source_log2fc']),source_insomnia_q=float(row['source_padj']),probe=probe,platform_symbol=platform_symbol,mapping=mapping,available_probes=len(options),frail_minus_nonfrail=float(beta[1]),hc3_se=se,p=p,same_sign=bool(float(row['source_log2fc'])*beta[1]>0)))
for row,q in zip(out,multipletests([x['p'] for x in out],method='fdr_bh')[1]):row['bh_q_panel']=float(q)
path=Path('outputs/stage6_spektor_gse40366_alias_sensitivity.csv')
with path.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
summary={'hgnc_sha256':__import__('hashlib').sha256(Path('work/geo_vitality/hgnc_complete_set.txt').read_bytes()).hexdigest(),'rescued_previous_symbol':[{k:x[k] for k in ('gene','platform_symbol','probe','p','bh_q_panel')} for x in out if x['mapping']=='previous_symbol'],'n_source_genes':len(source),'n_mapped':len(out),'n_unmapped':len(source)-len(out),'n_q_lt_05':sum(x['bh_q_panel']<.05 for x in out),'n_nominal_p_lt_05':sum(x['p']<.05 for x in out),'n_same_sign_all':sum(x['same_sign'] for x in out),'n_same_sign_q_lt_05':sum(x['same_sign'] and x['bh_q_panel']<.05 for x in out),'minimum_q':min(x['bh_q_panel'] for x in out),'n_frail':29,'n_nonfrail':38,'n_prefrail_excluded':84,'unmapped_genes':sorted(targets-{x['gene'] for x in out}),'source_platform':u,'method':'Exact GPL15988 symbol, or unambiguous HGNC previous symbol only when no exact probe exists; highest mean probe across all 151 categories; frail versus non-frail adjusted for sex with HC3 standard errors; BH across mapped panel','limitations':'Panel post-selected in another study, not preregistered; different populations, tissue and phenotype; bulk composition and unmeasured confounding; neither genetic nor score validation.'}
Path('outputs/stage6_spektor_gse40366_alias_sensitivity_receipt.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

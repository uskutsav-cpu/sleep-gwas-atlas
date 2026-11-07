import hashlib,json,math,re
from pathlib import Path

p=Path(__file__).resolve().parent; parent=p.parent
old=json.loads((parent/'prepare_receipt.json').read_text()); new=json.loads((p/'prepare_receipt.json').read_text()); n=new['selected_rows']
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
for k,path in old['output_paths'].items():assert sha(path)==old['output_sha256'][k]
for k,h in new['output_sha256'].items():assert sha(Path(new['work_path'])/f'{k}.sumstats.gz')==h
values={};deletes={};points={};corr=[]
for outcome in ['FI','FFS','HFRS','walking','grip']:
 for sleep in ['questionnaire','registry']:
  name=f'{outcome}_{sleep}'; folder=p if outcome in ('FI','FFS','HFRS') else parent
  log=folder/f'{name}.log';text=log.read_text();assert 'Analysis finished at' in text
  assert int(re.search(r'(\d+) SNPs with valid alleles',text).group(1))==n
  rg,se=map(float,re.search(r'Genetic Correlation: ([\d.-]+) \(([\d.]+)\)',text).groups())
  pv=float(re.search(r'Z-score: [\d.-]+\s+P: ([\d.eE-]+)',text).group(1))
  base=folder/f'{name}{outcome}.sumstats.gz_{sleep}.sumstats.gz'
  v={k:[float(x) for x in Path(str(base)+f'.{k}.delete').read_text().split()] for k in ['hsq1','hsq2','gencov']}
  assert all(len(x)==200 for x in v.values())
  assert all(math.isfinite(x) and x>0 for k in ('hsq1','hsq2') for x in v[k])
  d=[c/math.sqrt(a*b) for a,b,c in zip(v['hsq1'],v['hsq2'],v['gencov'])]; mean=sum(d)/200
  jse=math.sqrt(199/200*sum((x-mean)**2 for x in d));assert abs(jse-se)<0.00005
  values[name]=v;deletes[name]=d;points[name]=rg
  corr.append({'outcome':outcome,'insomnia':sleep,'rg':rg,'se':se,'p':pv,'bonferroni_10':min(1,10*pv),'log_sha256':sha(log)})
for sleep in ['questionnaire','registry']:
 assert all(values[f'{t}_{sleep}']['hsq2']==values[f'FI_{sleep}']['hsq2'] for t in ['FFS','HFRS','walking','grip'])
contrasts=[]
for t in ['FI','FFS','HFRS','walking','grip']:
 q=f'{t}_questionnaire';r=f'{t}_registry';assert values[q]['hsq1']==values[r]['hsq1']
 ds=[a-b for a,b in zip(deletes[q],deletes[r])];mean=sum(ds)/200
 se=math.sqrt(199/200*sum((x-mean)**2 for x in ds));diff=points[q]-points[r]
 pv=math.erfc(abs(diff/se)/math.sqrt(2))
 contrasts.append({'outcome':t,'questionnaire_minus_registry':diff,'paired_se':se,'unadjusted_95ci':[diff-1.95996398454*se,diff+1.95996398454*se],'p':pv,'bonferroni_5':min(1,5*pv)})
res={'marker_count':n,'correlations':corr,'contrasts':contrasts,'blocks':200,
     'matched_heritability_delete_vectors':True,'parent_receipt_sha256':sha(parent/'prepare_receipt.json'),
     'frailty_receipt_sha256':sha(p/'prepare_receipt.json'),'note':'Rounded point estimates and post-hoc families; upstream data and model assumptions remain.'}
(p/'results.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))

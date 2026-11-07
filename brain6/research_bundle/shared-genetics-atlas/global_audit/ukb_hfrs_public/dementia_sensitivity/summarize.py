import hashlib,json,math,re
from pathlib import Path

p=Path(__file__).resolve().parent; parent=p.parent
prep=json.loads((parent/'prepare_receipt.json').read_text()); stream=json.loads((p/'stream_receipt.json').read_text())
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
assert sha(stream['output_path'])==stream['output_sha256']
assert sha('/Volumes/Brain6Work/ukb_hfrs_public_20261028.sumstats.gz')==prep['ukb_hfrs_sha256']
for x in prep['sources'].values():assert sha(x['output_path'])==x['output_sha256']
assert stream['counts']['selected']==prep['selected_markers']
vals={};dels={};pts={};correlations=[]
for sleep in ('questionnaire','registry'):
 for kind in ('original','without'):
  folder=parent if kind=='original' else p
  name=f'ukb_{sleep}' if kind=='original' else f'without_{sleep}'
  prefix='ukb_hfrs_public_20261028' if kind=='original' else 'brain6-ukb-hfrs-without-dementia-20261028'
  log=folder/f'{name}.log';text=log.read_text();assert 'Analysis finished at' in text
  assert int(re.search(r'(\d+) SNPs with valid alleles',text).group(1))==prep['selected_markers']
  rg,se=map(float,re.search(r'Genetic Correlation: ([\d.-]+) \(([\d.]+)\)',text).groups())
  pv=float(re.search(r'Z-score: [\d.-]+\s+P: ([\d.eE-]+)',text).group(1))
  base=folder/f'{name}{prefix}.sumstats.gz_{sleep}.sumstats.gz'
  v={k:[float(x) for x in Path(str(base)+f'.{k}.delete').read_text().split()] for k in ('hsq1','hsq2','gencov')}
  assert all(len(x)==200 for x in v.values())
  assert all(math.isfinite(x) and x>0 for k in ('hsq1','hsq2') for x in v[k])
  ds=[c/math.sqrt(a*b) for a,b,c in zip(v['hsq1'],v['hsq2'],v['gencov'])];avg=sum(ds)/200
  secheck=math.sqrt(199/200*sum((x-avg)**2 for x in ds));assert abs(secheck-se)<0.00005
  key=f'{kind}_{sleep}';vals[key]=v;dels[key]=ds;pts[key]=rg
  correlations.append({'score':kind,'sleep':sleep,'rg':rg,'se':se,'p':pv,'bonferroni_4':min(1,4*pv),'log_sha256':sha(log)})
contrasts=[]
for sleep in ('questionnaire','registry'):
 a=f'original_{sleep}';b=f'without_{sleep}';assert vals[a]['hsq2']==vals[b]['hsq2']
 ds=[x-y for x,y in zip(dels[a],dels[b])];avg=sum(ds)/200;se=math.sqrt(199/200*sum((x-avg)**2 for x in ds))
 diff=pts[a]-pts[b];pv=math.erfc(abs(diff/se)/math.sqrt(2))
 contrasts.append({'sleep':sleep,'original_minus_without':diff,'paired_se':se,'p':pv,'bonferroni_2':min(1,2*pv),'ci95_unadjusted':[diff-1.95996398454*se,diff+1.95996398454*se]})
res={'markers':prep['selected_markers'],'correlations':correlations,'contrasts':contrasts,'blocks':200,'shared_insomnia_delete_vectors_match':True,'modified_source_sha256':stream['source_sha256'],'parent_prepare_sha256':sha(parent/'prepare_receipt.json'),'note':'Same cohort; rounded LDSC point estimates; per-file N discrepancy and overlap remain unresolved.'}
(p/'results.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))

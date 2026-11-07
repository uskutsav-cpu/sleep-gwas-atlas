import hashlib,json,math,re
from pathlib import Path

p=Path(__file__).resolve().parent; prep=json.loads((p/'prepare_receipt.json').read_text());n=prep['selected_markers']
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
assert sha('/Volumes/Brain6Work/ukb_hfrs_public_20261028.sumstats.gz')==prep['ukb_hfrs_sha256']
for x in prep['sources'].values():assert sha(x['output_path'])==x['output_sha256']
values={};dels={};pts={};correlations=[]
for h in ('ukb','finngen'):
 for sleep in ('questionnaire','registry'):
  name=f'{h}_{sleep}';log=p/f'{name}.log';text=log.read_text();assert 'Analysis finished at' in text
  assert int(re.search(r'(\d+) SNPs with valid alleles',text).group(1))==n
  rg,se=map(float,re.search(r'Genetic Correlation: ([\d.-]+) \(([\d.]+)\)',text).groups())
  pv=float(re.search(r'Z-score: [\d.-]+\s+P: ([\d.eE-]+)',text).group(1))
  hname='ukb_hfrs_public_20261028' if h=='ukb' else 'finngen_hfrs'
  base=p/f'{name}{hname}.sumstats.gz_{sleep}.sumstats.gz'
  v={k:[float(x) for x in Path(str(base)+f'.{k}.delete').read_text().split()] for k in ('hsq1','hsq2','gencov')}
  assert all(len(vv)==200 for vv in v.values())
  assert all(math.isfinite(x) and x>0 for k in ('hsq1','hsq2') for x in v[k])
  ds=[c/math.sqrt(a*b) for a,b,c in zip(v['hsq1'],v['hsq2'],v['gencov'])];avg=sum(ds)/200
  serr=math.sqrt(199/200*sum((x-avg)**2 for x in ds));assert abs(serr-se)<0.00005
  values[name]=v;dels[name]=ds;pts[name]=rg
  correlations.append({'hfrs':h,'sleep':sleep,'rg':rg,'se':se,'p':pv,'p_bonferroni_4':min(1,pv*4),'log_sha256':sha(log)})
for h in ('ukb','finngen'):assert values[f'{h}_questionnaire']['hsq1']==values[f'{h}_registry']['hsq1']
for sleep in ('questionnaire','registry'):assert values[f'ukb_{sleep}']['hsq2']==values[f'finngen_{sleep}']['hsq2']
spec=[('ukb_questionnaire','finngen_questionnaire'),('ukb_registry','finngen_registry'),('ukb_questionnaire','ukb_registry'),('finngen_questionnaire','finngen_registry')]
contrasts=[]
for a,b in spec:
 ds=[x-y for x,y in zip(dels[a],dels[b])];avg=sum(ds)/200;se=math.sqrt(199/200*sum((x-avg)**2 for x in ds))
 diff=pts[a]-pts[b];pv=math.erfc(abs(diff/se)/math.sqrt(2))
 contrasts.append({'contrast':f'{a} minus {b}','difference':diff,'paired_se':se,'p':pv,'p_bonferroni_4':min(1,4*pv),'ci95_unadjusted':[diff-1.95996398454*se,diff+1.95996398454*se]})
res={'markers':n,'correlations':correlations,'contrasts':contrasts,'blocks':200,'shared_trait_heritability_deletes_match':True,
     'prepare_receipt_sha256':sha(p/'prepare_receipt.json'),'stream_receipt_sha256':sha(p/'stream_receipt.json'),
     'note':'Post-hoc families; rounded LDSC points; cohort, phenotype, N discrepancy and reference limitations remain.'}
(p/'results.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))

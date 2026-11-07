from pathlib import Path
import csv,hashlib,json,re
p=Path(__file__).resolve().parent
atlas=next(q for q in p.parents if q.name=='shared-genetics-atlas')
preds=p.parent/'predictions.tsv'
root=atlas/'replication'/'mvp_frailty_genomewide'
log=root/'ldsc-20261029T052003-99497.log'
receipt=root/'stream-20261029T052003-99497.json'
prior=atlas/'replication'/'mvp_finngen_genomewide'/'stream_receipt.json'
assert json.loads(receipt.read_text())==json.loads(prior.read_text())
s=log.read_text(); assert 'Analysis finished' in s and s.count('1103239 SNPs with valid alleles.')==3
assert [s.index('/hfrs_ukb_sensitivity/'+x+'.sumstats.gz ...') for x in ('FI','FFS','HFRS')] == sorted(s.index('/hfrs_ukb_sensitivity/'+x+'.sumstats.gz ...') for x in ('FI','FFS','HFRS'))
matches=re.findall(r'Genetic Correlation: ([\d.\-]+) \(([\d.]+)\)\nZ-score: [\d.\-]+\nP: ([\deE.\-]+)',s)
assert len(matches)==3
outcomes=dict(zip(('FI','FFS','HFRS'),matches))
with preds.open(newline='') as f: rows=[r for r in csv.DictReader(f,delimiter='\t') if r['sleep_trait']=='insomnia']
assert len(rows)==3 and {r['target'] for r in rows}==set(outcomes)
models=('training_mean','knn3','ridge','bounded_ridge')
output=[]
for row in rows:
 target=row['target']; rg,se,pval=map(float,outcomes[target]); record={'target':target,'mvp_rg':rg,'mvp_se':se,'mvp_p':pval,'questionnaire_rg':float(row['observed'])}
 for model in models:
  record[model]=float(row[model]);record[model+'_absolute_error']=abs(float(row[model])-rg)
 output.append(record)
with (p/'comparison.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=output[0],delimiter='\t');w.writeheader();w.writerows(output)
summary={m:sum(r[m+'_absolute_error'] for r in output)/len(output) for m in models}
meta={'source_sha256':{str(x.relative_to(atlas)):hashlib.sha256(x.read_bytes()).hexdigest() for x in (preds,log,receipt,prior)},'mean_absolute_error':summary,'validation':'Outcome-substitution stress only; questionnaire features reused; no new fit.'}
(p/'receipt.json').write_text(json.dumps(meta,indent=2)+'\n')
print(json.dumps(output,indent=2));print('mean absolute errors:',summary)

import csv,json,math,pathlib,hashlib
p=pathlib.Path(__file__).resolve().parent;held=p.parent/'predictions.tsv';atlas=p.parents[3];src=atlas/'frailty_sensitivity/mak2023_cross_study'
with held.open() as f:r={x['target']:x for x in csv.DictReader(f,delimiter='\t') if x['held_sleep']=='insomnia'}
assert set(r)=={'FI','FFS','HFRS'}
rows=[];inputs=[held]
for source in ('ukb','r13','mvp'):
 for output,typ,old in [('FI','mak_FI','old_FI'),('FFS','mak_FP','old_FFS')]:
  f=src/f'mak_{source}_{typ}.exact.json';g=src/f'mak_{source}_{old}.exact.json';inputs.extend([f,g]);d=json.loads(f.read_text())[0];e=json.loads(g.read_text())[0];assert d['row_order_SNP_sha256']==e['row_order_SNP_sha256'] and d['n_snps']==e['n_snps'];x=r[output];target=float(d['rg'])
  row={'sleep_source':source,'target':typ,'reference_target':old,'n_snps':d['n_snps'],'target_rg':target,'target_se':d['rg_se'],'old_pair_rg':e['rg'],'original_model_target':float(x['observed'])}
  for model,key in [('ridge','ridge_prediction'),('mean','baseline'),('knn3','knn3')]:
   value=float(x[key]);row[model+'_prediction']=value;row[model+'_absolute_error']=abs(value-target)
  rows.append(row)
with (p/'transfer_comparison.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
summary=[]
for source in ('ukb','r13','mvp'):
 q=[x for x in rows if x['sleep_source']==source];assert len(q)==2
 summary.append({'sleep_source':source,'n_outcomes':2,**{m+'_mae':sum(x[m+'_absolute_error'] for x in q)/2 for m in ('ridge','mean','knn3')}})
with (p/'source_summary.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter='\t');w.writeheader();w.writerows(summary)
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
(p/'receipt.json').write_text(json.dumps({'input_sha256':{str(f):sha(f) for f in inputs},'rows':6,'prediction_reuse':'same UKB-insomnia-row-held-out predictions used for three source labels','limitations':'post-hoc two outcome proxy stress per source; reused/overlapping cohorts; not independent validation or causal prediction'},indent=2)+'\n');print(json.dumps(summary,indent=2))

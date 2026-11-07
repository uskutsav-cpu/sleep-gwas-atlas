import csv,pathlib,json,collections,math,statistics,hashlib
p=pathlib.Path('outputs/shared-genetics-atlas');o=p/'multiscale_embedding';models=o/'multistart_predictions.tsv';src=o/'source_matrix.tsv';
with models.open() as f:rows=[r for r in csv.DictReader(f,delimiter='\t') if r['sleep_trait']=='insomnia' and r['outcome'] in ('FI','FFS')]
with src.open() as f:source={(r['sleep_trait'],r['outcome']):r for r in csv.DictReader(f,delimiter='\t')}
assert len(rows)==20;refs={};paths={}
for name in ('FI','FFS'):
 q=p/'frailty_sensitivity'/(name+'_R13_ukbcommon.exact.json');paths[name]=q;a=json.loads(q.read_text());a=a[0] if isinstance(a,list) else a;assert a['n_snps']==1134953;refs[name]=a
out=[];methods=('column_mean','original_rank2','additive','multistart_rank2')
for r in rows:
 name=r['outcome'];ukb=float(source[('insomnia',name)]['rg']);assert float(r['observed'])==ukb
 for method in methods:
  val=float(r[method]);out.append({'repeat':r['repeat'],'fold':r['fold'],'outcome':name,'method':method,'UKB_insomnia_target':ukb,'R13_insomnia_target':refs[name]['rg'],'R13_target_se':refs[name]['rg_se'],'prediction':val,'error_vs_UKB':val-ukb,'error_vs_R13':val-refs[name]['rg']})
with (o/'registry_transfer_predictions.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
summary=[]
for name in ('FI','FFS'):
 for method in methods:
  a=[r for r in out if r['outcome']==name and r['method']==method];assert len(a)==10 and len({r['repeat'] for r in a})==10
  summary.append({'outcome':name,'method':method,'n_dependent_models':len(a),'UKB_insomnia_target':a[0]['UKB_insomnia_target'],'R13_insomnia_target':a[0]['R13_insomnia_target'],'R13_se':a[0]['R13_target_se'],'mean_prediction':statistics.mean(r['prediction'] for r in a),'prediction_min':min(r['prediction'] for r in a),'prediction_max':max(r['prediction'] for r in a),'RMSE_vs_UKB':math.sqrt(statistics.mean(r['error_vs_UKB']**2 for r in a)),'RMSE_vs_R13':math.sqrt(statistics.mean(r['error_vs_R13']**2 for r in a))})
with (o/'registry_transfer_summary.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter='\t');w.writeheader();w.writerows(summary)
(o/'registry_transfer_receipt.json').write_text(json.dumps({'source_sha256':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in (models,src,*paths.values())},'n_models_per_outcome':10,'interpretation':'post-hoc R13 sleep source shift with reused FI/FFS; not independent two-trait replication'},indent=2)+'\n')
for r in summary:print(r)

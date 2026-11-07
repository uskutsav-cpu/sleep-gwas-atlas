import csv,pathlib,collections,json,math,statistics,hashlib
root=pathlib.Path('outputs/shared-genetics-atlas/multiscale_embedding')
def read(name):
 with (root/name).open() as f:return list(csv.DictReader(f,delimiter='\t'))
src=read('source_matrix.tsv');pred=read('masked_predictions.tsv');lookup={(x['sleep_trait'],x['outcome']):x for x in src};assert len(src)==len(lookup)==96 and len(pred)==960
seen=collections.Counter();seen_pred=set()
for x in pred:
 k=(x['sleep_trait'],x['outcome']);assert k in lookup and float(x['observed'])==float(lookup[k]['rg']);rep=(x['repeat'],*k);assert rep not in seen_pred;seen_pred.add(rep);seen[k]+=1;x['se']=float(lookup[k]['se']);assert x['se']>0
assert set(seen.values())=={10}
groups={'all':pred,'brain': [x for x in pred if x['outcome'] not in ('FI','FFS')], 'FI':[x for x in pred if x['outcome']=='FI'],'FFS':[x for x in pred if x['outcome']=='FFS'],'exclude_longsleep':[x for x in pred if x['sleep_trait']!='longsleep']}
out=[]
for g,rows in groups.items():
 for name in ('column_mean','rank2','rank2_clipped','additive'):
  errors=[];scaled=[]
  for x in rows:
   p=float(x['rank2'] if name=='rank2_clipped' else x[name]);p=min(1.,max(-1.,p)) if name=='rank2_clipped' else p;e=p-float(x['observed']);errors.append(e);scaled.append(e/x['se'])
  out.append({'group':g,'method':name,'n_predictions_dependent':len(rows),'n_unique_cells':len({(x['sleep_trait'],x['outcome']) for x in rows}),'rmse':math.sqrt(sum(e*e for e in errors)/len(errors)),'median_absolute_error':statistics.median(abs(e) for e in errors),'se_scaled_rmse_descriptive':math.sqrt(sum(e*e for e in scaled)/len(scaled))})
with (root/'error_audit.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
meta={'source_hashes':{n:hashlib.sha256((root/n).read_bytes()).hexdigest() for n in ('source_matrix.tsv','masked_predictions.tsv')},'n_predictions':len(pred),'n_unique_cells':len(lookup),'rank2_out_of_range':sum(abs(float(x['rank2']))>1 for x in pred),'caveat':'post hoc repeated dependent predictions; SE scaling ignores prediction variance and covariance'}
(root/'error_audit_receipt.json').write_text(json.dumps(meta,indent=2)+'\n')
for x in out:print(x)

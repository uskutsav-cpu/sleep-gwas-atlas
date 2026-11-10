import csv,json,math,hashlib
from pathlib import Path

def read(p):
 with open(p,newline='') as f:return list(csv.DictReader(f))
target=read('outputs/stage6_kunorozva_meta_candidate_extract.csv');t={r['candidate']:r for r in target}
assert len(t)==len(target)==33 and all(r['candidate_pair_match']=='True' for r in target)
configs=[('15_weighted','outputs/unvalidated_15_variant_score_prototype.csv',[('FI','FI_beta_weight'),('FFS','FFS_beta_weight_inferred'),('insomnia','insomnia_beta_weight')]),('22_weighted','outputs/unvalidated_22_variant_score_prototype.csv',[('FI','FI_beta_weight'),('FFS','FFS_beta_weight_inferred'),('insomnia','UKB_insomnia_beta_weight')]),('22_count','outputs/unvalidated_22_shared_direction_score.csv',[('count','weight')])]
out=[]
for score,path,weights in configs:
 rows=read(path)
 for x in rows:
  r=t[x['variant']];assert {x['counted_allele'],x['other_allele']}=={r['ref'],r['alt']}
  # Conditional on source beta=ALT effect and a valid biallelic diploid contrast.
  aligned=float(r['beta'])*(1 if x['counted_allele']==r['alt'] else -1)
  for label,col in weights:
   wt=float(x[col]);assert wt!=0 and aligned!=0
   out.append({'score':score,'weight':label,'variant':x['variant'],'counted_allele':x['counted_allele'],'source_alt':r['alt'],'source_beta_alt_assumed':float(r['beta']),'counted_allele_beta_conditional':aligned,'source_se':float(r['stderr_beta']),'source_neg_log10p':float(r['neg_log_pvalue']),'weight_value':wt,'direction_matches':int(wt*aligned>0)})
p=Path('outputs/stage6_kunorozva_meta_direction_comparison.csv')
with p.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
summary={}
for score,path,weights in configs:
 for label,_ in weights:
  x=[r for r in out if r['score']==score and r['weight']==label]
  summary[f'{score}_{label}']={'n':len(x),'direction_matches':sum(r['direction_matches'] for r in x),'nominal_p_lt_05':sum(r['source_neg_log10p']>-math.log10(.05) for r in x)}
receipt={'summary':summary,'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'limitations':'Postselected, same or overlapping discovery cohorts, phenotype differences and conditional ALT effect and biallelic dosage assumptions; no independent validation or predictions.'}
Path('outputs/stage6_kunorozva_meta_direction_comparison_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

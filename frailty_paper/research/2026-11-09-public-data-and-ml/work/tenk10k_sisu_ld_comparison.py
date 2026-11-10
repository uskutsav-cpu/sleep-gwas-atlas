import requests,json,hashlib,numpy as np,pandas as pd
from pathlib import Path
out=Path('outputs'); anchors=['5:136187682:T:G','5:136235970:T:C','5:136131620:C:G']; ref=pd.read_csv(out/'stage6_tenk10k_expanded_eur_ld_order.csv'); d=pd.read_csv('work/stage6_tenk10k_susie_expanded/aligned.csv'); z=np.load(out/'stage6_tenk10k_expanded_eur_ld.npz'); matrix=z['r']; idx=dict(zip(ref.variant_GRCh38,range(len(ref)))); mask=d.variant_GRCh38.to_list(); assert all(x in idx for x in anchors+mask); summaries=[]
for k,a in enumerate(anchors):
 p={'variant':a,'window':1000000,'panel':'sisu42','r2_thresh':0};resp=requests.get('https://api.finngen.fi/api/ld',params=p,timeout=120);resp.raise_for_status();j=resp.json(); assert all(x['variation1']==a for x in j['ld']); path=out/f'stage6_tenk10k_sisu_anchor_{k+1}.json';path.write_text(json.dumps(j,separators=(',',':'))); vals={x['variation2']:x['r2'] for x in j['ld']};rows=[]
 for v in mask:
  if v in vals:
   e=float(matrix[idx[a],idx[v]])**2;s=float(vals[v]);rows.append({'anchor':a,'variant':v,'sisu_r2':s,'eur_r2':e,'abs_difference':abs(e-s)})
 df=pd.DataFrame(rows); df.to_csv(out/f'stage6_tenk10k_sisu_anchor_{k+1}_comparison.csv.gz',index=False)
 print(a,len(j['ld']),len(df),'targets',[(v,vals.get(v)) for v in anchors],flush=True)
 summaries.append({'anchor':a,'url':resp.url,'response_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'response_rows':len(j['ld']),'matched_fi_snp_count':len(df),'fi_snp_count':len(mask),'median_absolute_difference':float(df.abs_difference.median()),'p95_absolute_difference':float(df.abs_difference.quantile(.95)),'count_difference_gt_0_1':int((df.abs_difference>.1).sum()),'targets':{v:vals.get(v) for v in anchors}})
(out/'stage6_tenk10k_sisu_ld_comparison.json').write_text(json.dumps(summaries,indent=2)+'\n')

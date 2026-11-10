import json,csv,glob,zipfile,datetime
from pathlib import Path
files=glob.glob('work/stage5_spliceai/*.json'); rows=[]
for f in files:
 d=json.load(open(f)); token=Path(f).name.split('_')[1].split('.')[0]; variant=f"{d.get('chrom',d['variant'].split('-')[0])}:{d.get('pos',d['variant'].split('-')[1])}:{d.get('ref',d['variant'].split('-')[2])}>{d.get('alt',d['variant'].split('-')[3])}"
 variant = variant if variant.startswith('chr') else 'chr'+variant
 base={'variant_GRCh38':variant,'model':token,'source':d.get('source'),'window_bp':d.get('distance'),'mask':d.get('mask'),'error':d.get('error',''),'retrieved_utc':'2026-10-27','api_url':f"https://{token}-38-xwkwwwxdwq-uc.a.run.app/{token}/?hg=38&variant={d['variant']}"}
 scores=d.get('scores') or []
 if not scores: rows.append(base)
 for s in scores:
  r=base.copy()
  for k in ['g_id','g_name','t_id','t_priority','t_strand','t_type','DS_AG','DS_AL','DS_DG','DS_DL','DP_AG','DP_AL','DP_DG','DP_DL','DS_SG','DS_SL','DP_SG','DP_SL']:r[k]=s.get(k,'')
  rows.append(r)
fields=['variant_GRCh38','model','source','window_bp','mask','g_id','g_name','t_id','t_priority','t_strand','t_type','DS_AG','DS_AL','DS_DG','DS_DL','DP_AG','DP_AL','DP_DG','DP_DL','DS_SG','DS_SL','DP_SG','DP_SL','error','retrieved_utc','api_url']
with open('outputs/stage5_priority_splice_model_predictions.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
with zipfile.ZipFile('outputs/stage5_priority_splice_model_raw_responses.zip','w',zipfile.ZIP_DEFLATED) as z:
 for f in files:z.write(f,Path(f).name)
print(len(files),'responses,',len(rows),'rows')

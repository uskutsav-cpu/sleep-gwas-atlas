#!/usr/bin/env python3
"""Retrieve browser-visible coloc pairs for fixed ±250 kb marker windows."""
import csv,json,pathlib,requests,hashlib,zipfile,time
from urllib.parse import quote
with open('outputs/stage5_score_queue.csv') as f:queue=[x for x in csv.DictReader(f) if 'exploratory_marker' in x['input_classes']]
assert len(queue)==23
records=[];receipt=[];source=pathlib.Path('work/stage5_cascade_coloc_raw');source.mkdir(exist_ok=True)
for x in queue:
 pos=int(x['position_GRCh38']);start=max(0,pos-250000);end=pos+250000;region=f"{x['chr']}:{start}-{end}";url='https://cascade.finngen.fi/api/cascade/coloc/by_region/'+quote(region,safe='')
 for attempt in range(4):
  try:
   r=requests.get(url,timeout=45);r.raise_for_status();d=r.json()
   if 'pairs' not in d:raise ValueError(str(d.get('reason','missing pairs')))
   break
  except Exception:
   if attempt==3:raise
   time.sleep(2**attempt)
 raw=r.content;name=hashlib.sha256(x['variant_key'].encode()).hexdigest()[:20]+'.json';(source/name).write_bytes(raw);receipt.append({'variant':x['variant_key'],'region':region,'url':url,'file':name,'sha256':hashlib.sha256(raw).hexdigest(),'pairs':len(d['pairs'])})
 for p in d['pairs']:records.append({'index_marker':x['variant_key'],'marker_labels':x['input_labels'],'queried_region':region,**p})
 print(x['variant_key'],'pairs',len(d['pairs']),flush=True)
cols=sorted({k for r in records for k in r});
with open('outputs/stage5_cascade_regional_coloc_pairs.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(records)
with zipfile.ZipFile('outputs/stage5_cascade_regional_coloc_raw.zip','w',zipfile.ZIP_DEFLATED) as z:
 for x in receipt:z.write(source/x['file'],x['file'])
pathlib.Path('outputs/stage5_cascade_regional_coloc_receipt.json').write_text(json.dumps({'source':'public CASCADE by-region API','windows':'23 markers, ±250 kb, zero-based endpoint convention unverified','note':'Browser displays selected colocalization results; no entry does not prove no colocalization. Windows may overlap and repeat same pairs.','requests':receipt},indent=2)+'\n')
print('total rows',len(records),'markers',len(queue))

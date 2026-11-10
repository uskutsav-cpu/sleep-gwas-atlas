#!/usr/bin/env python3
from pathlib import Path
import requests,json,hashlib
root=Path(__file__).resolve().parents[1];o=root/'outputs';anchors={'rs35624276':'7:114400048:G:A','rs12536335':'7:114403104:A:G','rs62477310':'7:115311487:T:C','rs2396682':'7:113765895:T:G'};receipt=[]
for rs,v in anchors.items():
 u='https://api.finngen.fi/api/ld';params={'variant':v,'window':2000000,'panel':'sisu42','r2_thresh':0};r=requests.get(u,params=params,timeout=60);r.raise_for_status();j=r.json();assert isinstance(j.get('ld'),list) and all(x['variation1']==v for x in j['ld']);p=o/f'stage6_foxp2_sisu_ld_{rs}.json';p.write_text(json.dumps(j,separators=(',',':')));receipt.append({'rsid':rs,'variant_grch38':v,'url':r.url,'n_rows':len(j['ld']),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()});print(rs,len(j['ld']),flush=True)
(o/'stage6_foxp2_sisu_anchor_receipt.json').write_text(json.dumps({'panel':'sisu42','window_bp':2000000,'threshold':0,'queries':receipt,'limitations':'API r2 and D-prime are not a signed LD correlation matrix; reference panel is not in-sample GWAS LD.'},indent=2)+'\n')

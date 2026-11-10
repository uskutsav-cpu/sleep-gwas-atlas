#!/usr/bin/env python3
"""Extract exact-position candidate rows from the official BBJ insomnia autosomal TSV.
Input: gzip-decompressed TSV on stdin; outputs are written only after EOF.
"""
import csv, sys, json
from pathlib import Path
p=Path(__file__).resolve().parent
with (p/'stage6_japanese_knee_candidate_matches.csv').open(newline='') as f:
 candidates=list(csv.DictReader(f))
lookup={}
for x in candidates:
 lookup.setdefault((x['chromosome'],x['position']),[]).append(x)
r=csv.DictReader(sys.stdin,delimiter='\t'); rows=[]; count=0
for x in r:
 count+=1
 for c in lookup.get((x['CHR'],x['POS']),[]):
  exact={x['Allele1'].upper(),x['Allele2'].upper()}=={c['counted_allele'].upper(),c['other_allele'].upper()}
  rows.append({'score':c['score'],'candidate_id':c['variant'],'counted_allele':c['counted_allele'],'other_allele':c['other_allele'],'exact_pair':int(exact),**x})
fields=['score','candidate_id','counted_allele','other_allele','exact_pair']+r.fieldnames
out=p/'stage6_bbj_insomnia_candidate_matches.csv'
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
(p/'stage6_bbj_insomnia_extract_receipt.json').write_text(json.dumps({'source':'https://humandbs.dbcls.jp/files/hum0197/hum0197.v3.BBJ.Ins.v1.zip','autosomal_stream_rows':count,'candidate_rows_at_position':len(rows),'exact_pair_rows':sum(x['exact_pair'] for x in rows),'note':'Completed decompressed stream; candidate source is stage6_japanese_knee_candidate_matches.csv. Verify gzip stream status separately.'},indent=2)+'\n')
print('rows scanned',count,'matches',len(rows),'exact pairs',sum(x['exact_pair'] for x in rows),flush=True)

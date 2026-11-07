#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
with (p/'exact_qtl_rows.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
r=json.loads((p/'receipt.json').read_text())
assert len(rows)==2 and {x['rsid'] for x in rows}=={'rs6874138','rs325485'}
assert all(x['context']=='ast' and x['target_ensembl_id']=='ENSG00000112874' and float(x['pval'])>.0125 and 'within_cis_window=false' in x['target_info'] for x in rows)
for ctx in ('ast','opc'):
 for rs in ('rs325485','rs6874138'):
  assert r['slots'][ctx+'_'+rs]['row_count']==sum(x['context']==ctx and x['rsid']==rs for x in rows)
for rs in ('rs325485','rs6874138'):
 data=(p/(rs+'_portal_response')).read_bytes();assert hashlib.sha256(data).hexdigest()==r['portal'][rs]['sha256'] and b'No records found for' in data and rs.encode() in data
assert hashlib.sha256((p/'exact_qtl_rows.tsv').read_bytes()).hexdigest()==r['rows_sha256']
receipt={'returned_rows':len(rows),'slots_reconciled':4,'portal_empty_responses_confirmed':2,'within_cis_window_false':2,'row_sha256':r['rows_sha256'],'limitations':'Reconciles extracted records and saved portal responses; does not prove that missing associations were tested, validate the source statistical model, or establish target mediation.'}
(p/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

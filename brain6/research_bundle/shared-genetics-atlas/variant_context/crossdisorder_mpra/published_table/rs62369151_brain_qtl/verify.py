#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
with (p/'exact_qtl_rows.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
with (p/'portal_response').open() as f:portal=list(csv.DictReader(f,delimiter='\t'))
r=json.loads((p/'receipt.json').read_text());assert len(rows)==4 and len({(x['context'],x['rsid']) for x in rows})==4
assert all(float(x['pval'])<.0125 and float(x['FDR'])<.05 for x in rows)
for x in rows:
 assert r['slots'][x['context']+'_'+x['rsid']]['row_count']==1
 if x['rsid']=='rs62369151':
  hits=[a for a in portal if a['variantId']==x['variant_id'] and a['target']==x['target_ensembl_id'] and a['context']==x['context'] and a['study']=='ROSMAP_CUIMC1_2_MIT' and a['xQTLtype']=='snuc-eQTL']
  assert len(hits)==1
  for source,target in [('beta','beta_non_ref'),('se','beta_se_non_ref'),('pvalue','pval'),('FDR','FDR')]:assert math.isclose(float(hits[0][source]),float(x[target]),rel_tol=1e-9,abs_tol=1e-18)
assert hashlib.sha256((p/'exact_qtl_rows.tsv').read_bytes()).hexdigest()==r['rows_sha256']
assert hashlib.sha256((p/'portal_response').read_bytes()).hexdigest()==r['portal']['sha256']
receipt={'four_slots_reconciled':True,'rs62369151_portal_rows_reproduced':2,'all_four_below_postselected_0125_threshold':True,'rows_sha256':r['rows_sha256'],'limitations':'Checks row matching and arithmetic, not independence, causal mechanisms, disease association colocalization, or statistical calibration.'}
(p/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

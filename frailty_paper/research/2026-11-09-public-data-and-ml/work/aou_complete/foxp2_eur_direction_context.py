import csv,json,hashlib
from pathlib import Path
sources=[Path('outputs/stage6_aou_eur_candidate_association.csv'),Path('outputs/stage6_aou_eur_walking_candidate_association.csv'),Path('outputs/stage6_aou_eur_fatigue_candidate_association.csv')]
wanted={'NS_333.2','NS_325.2','SS_807','SS_807.21'};rows=[]
for path in sources:
 with path.open(newline='') as f:
  for r in csv.DictReader(f):
   if r['candidate']=='rs12536335' and r['phenotype'] in wanted:
    assert r['status']=='matched' and r['ancestry']=='eur' and r['sequencing_type']=='genome' and (r['chrom'],r['position'],r['ref'],r['alt'])==('7','114403104','A','G')
    rows.append({k:r[k] for k in ('candidate','phenotype','chrom','position','ref','alt','beta','se','pvalue','af')})
assert len(rows)==4 and {r['phenotype'] for r in rows}==wanted
p=Path('outputs/stage6_foxp2_aou_eur_direction_context.csv')
with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
rec={'inputs_sha256':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in sources},'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'rows':len(rows),'note':'Marginal associations, different phenotypes; ALT effect coding is supported by a public-code audit but needs verification.'}
Path('outputs/stage6_foxp2_aou_eur_direction_context_receipt.json').write_text(json.dumps(rec,indent=2)+'\n');print('four exact rows extracted; sha256',rec['output_sha256'])

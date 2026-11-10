#!/usr/bin/env python3
import csv,hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parent.parent;out=root/'outputs'
spec=[('15','reference_15_score_variance_decomposition.csv','reference_counted_allele_frequency'),
      ('22','22_score_reference_group_allele_frequencies.csv','reference_EUR_frequency')]
report={}
for label,name,column in spec:
 path=out/name
 with path.open() as f:rows=list(csv.DictReader(f))
 assert len(rows)==int(label) and len({r['variant'] for r in rows})==len(rows)
 vals=sorted((min(float(r[column]),1-float(r[column])),r['variant']) for r in rows)
 report[label]={'n_markers':len(vals),'at_or_below_1pct':sum(x<=.01 for x,_ in vals),
                'minimum_maf':vals[0][0],'minimum_marker':vals[0][1],
                'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
receipt={'panel':'1000 Genomes Phase 3 503 European-labeled samples; previously computed reference frequencies',
         'results':report,'limitation':'External reference allele-frequency summary only; not NCDS genotype or imputation coverage. Multiallelic contrasts and target QC unverified.'}
(out/'stage6_ncds_reference_maf_check.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(report,indent=2))

#!/usr/bin/env python3
"""Screen the downloaded public CKB PheWeb phenotype-code inventory."""
import csv,hashlib,json,re
from pathlib import Path
root=Path(__file__).resolve().parent.parent;file=root/'work/stage6_ckb_phenotypes.tsv';out=root/'outputs'
with file.open(newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
pat=r'frail|grip|sleep|insomn|walk|weak';codes=['g47','f51','r54','r53','m62','m62.81','r26','r26.2'];positives=['e11','i10']
labels=[x['phenocode'] for x in rows];result={'source':'https://pheweb.ckbiobank.org/download/phenotypes.tsv','sha256':hashlib.sha256(file.read_bytes()).hexdigest(),'rows':len(rows),'unique_codes':len(set(x.lower() for x in labels)),'rows_with_case_count':sum(bool(x['num_cases']) for x in rows),'keyword_pattern':pat,'keyword_matches':sorted(x for x in labels if re.search(pat,x,re.I)),'exact_codes':{q:sum(x.lower()==q for x in labels) for q in codes},'positive_controls':{q:sum(x.lower()==q for x in labels) for q in positives},'limits':'Searches only listed phenocode values and selected exact codes. Codes can be abbreviated or differently mapped; missing code does not establish absence of a trait or an external GWAS, no individual validation or score association.'}
(out/'stage6_ckb_catalog_screen.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['rows','unique_codes','rows_with_case_count','keyword_matches','exact_codes','positive_controls']},indent=2))

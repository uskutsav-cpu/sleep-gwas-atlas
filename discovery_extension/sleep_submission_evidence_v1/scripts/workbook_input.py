#!/usr/bin/env python3
"""Typed interchange for editable TSV-derived supplementary workbook."""
import csv, json, re, sys
from pathlib import Path
P=Path(__file__).resolve().parents[1]
FILES=[('Replicated pairs','replicated_23.tsv'),('Global correlations','global_1200.tsv'),
       ('Replication family','replication_family_217.tsv'),('Heterogeneity','heterogeneity_7.tsv'),
       ('Sources','phenotype_and_source_metadata.tsv'),('Sleep definitions','sleep_trait_metadata.tsv'),
       ('Novelty','NOVELTY_MASTER_1200.tsv'),('Sensitivity','sensitivity_results.tsv')]
out=[]
for sheet,file in FILES:
    with (P/'tables'/file).open(newline='') as f:
        reader=csv.DictReader(f,delimiter='\t');heads=reader.fieldnames;rows=[]
        for record in reader:
            values=[]
            for h in heads:
                v=record[h]
                if v in ('True','False'): v=v=='True'
                elif not re.search('(^|_)(id|accession|PMID|pmid|DOI|doi|code)($|_)',h) and re.fullmatch(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?',v):
                    v=float(v) if any(c in v.lower() for c in '.e') else int(v)
                values.append(v)
            rows.append(values)
    out.append({'sheet':sheet,'source':'tables/'+file,'headers':heads,'rows':rows})
dest=Path(sys.argv[1]);dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(json.dumps(out))
print(json.dumps({'worksheets':len(out),'rows':sum(len(x['rows']) for x in out)}))

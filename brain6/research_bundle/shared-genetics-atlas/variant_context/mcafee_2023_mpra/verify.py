#!/usr/bin/env python3
import csv,gzip,hashlib,json
from pathlib import Path
p=Path(__file__).resolve().parent;r=json.loads((p/'lookup_receipt.json').read_text())
with gzip.open(p/'cre_sequences.tsv.gz','rt') as f: rows=list(csv.reader(f,delimiter='\t'))
assert rows[0]==['names','seq'] and len(rows)-1==10346
names=[x[0] for x in rows[1:]];assert len(set(names))==len(names)
refs={x[:-8] for x in names if x.endswith('_scz_ref')};alts={x[:-8] for x in names if x.endswith('_scz_alt')}
assert len(refs)==5173 and refs==alts
with (p.parents[1]/'prediction_rank_comparison/variant_ranks.tsv').open() as f:candidates={x['snp'] for x in csv.DictReader(f,delimiter='\t')}
with (p.parent/'rs62264764_dense_qtl/multisignal_pilot/variants.tsv').open() as f:region={x['rsid'] for x in csv.DictReader(f,delimiter='\t')}
assert not refs & (candidates|{'rs62264767'}) and not refs & region
assert r['candidate_matches']==r['regional_matches']=={}
receipt={'paired_variant_identifiers':len(refs),'constructs':len(names),'candidate_intersection_count':0,'regional_intersection_count':0,'lookup_receipt_sha256':hashlib.sha256((p/'lookup_receipt.json').read_bytes()).hexdigest(),'scope':'Independent name-level coverage check only; no test of genomic allele identity or biological absence.'}
(p/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

#!/usr/bin/env python3
import csv,hashlib,json
from pathlib import Path
p=Path(__file__).resolve().parent;root=p.parent.parent
prov=json.loads((root.parent/'brain6-alphagenome-pilot/results/provenance.json').read_text())['input_sha256']
src=Path('/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v2/inputs')
def read(x):
 with x.open() as f:return list(csv.DictReader(f,delimiter='\t'))
candidates=[x for x in read(root/'evidence/BRAIN6_FINAL_CANDIDATES.tsv') if x['pair_id']=='insomnia__adhd'];observed={x['variant']:x for x in read(p/'indexed_results.tsv')}
assert len(candidates)==6
for c in candidates:
 rows=[]
 for trait in ('insomnia','adhd'):
  name=c['candidate_locus_id']+'__'+trait+'.tsv';path=src/name
  assert hashlib.sha256(path.read_bytes()).hexdigest()==prov[name]
  matched=[x for x in read(path) if x['SNP']==c['lead_variants']];assert len(matched)==1;rows.append(matched[0])
 a,b=rows;r=observed[c['lead_variants']];assert a['A1']==r['effect_allele'] and {a['A1'],a['A2']}=={b['A1'],b['A2']}
 beta=float(b['BETA'])*(1 if b['A1']==a['A1'] else -1)
 status='not estimable' if r['status']!='MATCHED' else 'same' if beta*float(r['beta_aligned'])>0 else 'opposite or zero'
 print(c['lead_variants'],a['A1'],beta,status)

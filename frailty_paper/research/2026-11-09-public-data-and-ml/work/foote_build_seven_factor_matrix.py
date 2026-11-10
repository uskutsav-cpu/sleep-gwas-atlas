import csv,json,hashlib
from pathlib import Path
factors={46:'general',47:'limited_social_support',48:'unhealthy_lifestyle',49:'multimorbidity',50:'metabolic_problems',51:'poorer_cognition',52:'disability'}
source=list(csv.DictReader(open('outputs/shared22_eight_source_evidence_matrix.csv')))
results={}
for i,name in factors.items():
 acc=f'GCST906240{i}'; data=json.load(open(f'outputs/foote_{acc}_indexed_candidate_lookup.json')); assert data['accession']==acc
 if i in (47,48,50,51):
  idx=Path(f'work/{acc}.h.tsv.gz.tbi'); md=Path(f'work/{acc}.md5sum.txt').read_text(); exp=next(line.split()[0] for line in md.splitlines() if line.endswith('.tbi')); assert hashlib.md5(idx.read_bytes()).hexdigest()==exp==data['index_md5']
 results[name]={r['variant']:r for r in data['factor_candidate_rows']}
 assert len(results[name])==22
out=[]
for s in source:
 row={k:s[k] for k in ('variant','chr_GRCh38','position_GRCh38','score_counted_allele','other_allele','nearest_gene_body_positional')}
 for name in factors.values():
  r=results[name][s['variant']]; assert r['score_allele']==s['score_counted_allele'] and r['other_allele']==s['other_allele'] and r['position_GRCh38']==int(s['position_GRCh38']);
  for k in ('matched_rows','beta','se','p','q_metric_p'): row[f'{name}_{k}']=r.get(k,'')
 out.append(row)
with open('outputs/foote_seven_factors_shared22_lookup.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0].keys());w.writeheader();w.writerows(out)

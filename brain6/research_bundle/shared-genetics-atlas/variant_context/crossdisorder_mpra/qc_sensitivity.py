import csv,json,hashlib
from pathlib import Path
p=Path(__file__).resolve().parent;root=Path('/Volumes/Extreme SSD/brain6-crossdisorder-mpra')
names=['v124','v124_all_listed','v124_prop30']
def read(name):
 with (root/f'limma_pilot_{name}_all.tsv').open() as f:return {r['variant']:r for r in csv.DictReader(f,delimiter='\t')}
models={n:read(n) for n in names}
with (p/'ranked_coverage.tsv').open() as f:covered=[r for r in csv.DictReader(f,delimiter='\t') if r['dna_present']=='True' and r['rna_present']=='True']
rows=[]
for r in covered:
 item={k:r[k] for k in ['rsid','region','source_key','listed_excluded']}
 for n in names:
  v=models[n].get(r['source_key'])
  for col in ['logFC','P.Value','adj.P.Val']:item[n+'_'+col]=v[col] if v else 'EXCLUDED'
 rows.append(item)
assert len(rows)==42
with (p/'qc_sensitivity_42.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
common=[r for r in rows if all(r[n+'_adj.P.Val']!='EXCLUDED' for n in names)];assert len(common)==41
flips={n:sum((float(r['v124_adj.P.Val'])<.05)!=(float(r[n+'_adj.P.Val'])<.05) for r in common) for n in names[1:]}
sha=lambda x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
receipt={'models':{n:{'n':len(v),'fdr_lt_0.05':sum(float(r['adj.P.Val'])<.05 for r in v.values()),'sha256':sha(root/f'limma_pilot_{n}_all.tsv')} for n,v in models.items()},'matched_variants':len(rows),'common_variants':len(common),'threshold_changes_from_no_filter':flips,'output_sha256':sha(p/'qc_sensitivity_42.tsv'),'code_sha256':sha(p/'limma_pilot.R')}
(p/'qc_sensitivity_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

import csv,openpyxl,json
from pathlib import Path
with open('outputs/stage5_variant_gene_tissue_phenotype_integrated.csv') as f:g={r['gene_symbol'] for r in csv.DictReader(f)}
w=openpyxl.load_workbook('work/hfrs_2025_tables.xlsx',read_only=True,data_only=True);out={}
for n,col in [('ST11',10),('ST12',10),('ST13',0),('ST14',0)]:
 rows=[r for r in list(w[n].values)[2:] if len(r)>col and r[col] in g]
 out[n]=[list(r[:11] if col==10 else r[:7]) for r in rows]
Path('outputs/stage6_hfrs_candidate_pqtl_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print({k:len(v) for k,v in out.items()})

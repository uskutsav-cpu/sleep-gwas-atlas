import csv
from pathlib import Path
cand={}
for x in csv.DictReader(open('outputs/stage3_initial_variant_gene_table.csv')):
 g=x['gene ID'].split('.')[0]
 if g.startswith('ENSG'): cand.setdefault(g,{'markers':set(),'symbols':set()});cand[g]['markers'].add(x['Variant']);cand[g]['symbols'].add(x['Candidate gene'])
a={x['GENE']:x for x in csv.DictReader(open('outputs/stage4_mak_fi_magma_gene_results.csv'))}
b={x['GENE']:x for x in csv.DictReader(open('outputs/stage4_ye_ffs_magma_gene_results.csv'))}
rows=[]
for g,v in sorted(cand.items()):
 x=a.get(g,{});y=b.get(g,{});rows.append({'gene_id':g,'candidate_symbols':'|'.join(sorted(v['symbols'])),'candidate_markers':'|'.join(sorted(v['markers'])),'mak_fi_gene_p':x.get('P',''),'mak_fi_gene_BH_q_19129':x.get('bh_q_all_genes',''),'ye_ffs_gene_p':y.get('P',''),'ye_ffs_gene_BH_q_19015':y.get('bh_q_all_19015_genes',''),'interpretation':'Two separate LD-aware gene tests; no formal difference test, shared-variant or mediation inference'})
with open('outputs/stage4_candidate_standard_fi_ffs_gene_comparison.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print('unique candidates',len(rows),'both tested',sum(bool(r['mak_fi_gene_p']) and bool(r['ye_ffs_gene_p']) for r in rows),'both q<.05',sum(bool(r['mak_fi_gene_p']) and bool(r['ye_ffs_gene_p']) and float(r['mak_fi_gene_BH_q_19129'])<.05 and float(r['ye_ffs_gene_BH_q_19015'])<.05 for r in rows))
for r in rows:
 if r['mak_fi_gene_p'] and r['ye_ffs_gene_p'] and float(r['mak_fi_gene_BH_q_19129'])<.05 and float(r['ye_ffs_gene_BH_q_19015'])<.05:print(r['gene_id'],r['candidate_symbols'],r['candidate_markers'],r['mak_fi_gene_BH_q_19129'],r['ye_ffs_gene_BH_q_19015'])

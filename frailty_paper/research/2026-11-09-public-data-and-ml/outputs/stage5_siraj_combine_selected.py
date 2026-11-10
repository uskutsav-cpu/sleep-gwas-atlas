#!/usr/bin/env python3
"""Combine selected analysis rows, preserving repeated phenotype/gene records."""
import csv,argparse
p=argparse.ArgumentParser();p.add_argument('--gtex',required=True);p.add_argument('--traits',required=True);a=p.parse_args()
with open('outputs/stage5_siraj2026_neural_mpra_exact_alleles.csv') as f:m={f"chr{x['variant_GRCh37_source']}":x['variant_GRCh38'] for x in csv.DictReader(f)}
rows=[]
for path,source in [(a.gtex,'GTEx'),(a.traits,'complex_trait')]:
 with open(path) as f:
  for r in csv.DictReader(f):
   r['variant_GRCh38']=m[r['variant']];r['analysis_source']=source;rows.append(r)
fields=['variant_GRCh38','variant','variant_hg38','analysis_source','cell_type','library','best_library','type','gene','tissue','trait','cohort','pip','ID','A_log2FC','A_logPadj_BF','B_log2FC','B_logPadj_BF','log2FC','logPadj_BF','Log2Skew','Skew_logPadj','active_screen','emvar_screen']
with open('outputs/stage5_siraj2026_selected_library_records.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)
print('rows',len(rows),'variants',len({r['variant_GRCh38'] for r in rows}))

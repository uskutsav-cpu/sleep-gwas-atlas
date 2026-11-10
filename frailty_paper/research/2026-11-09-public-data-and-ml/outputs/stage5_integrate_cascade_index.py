#!/usr/bin/env python3
"""Add CASCADE coverage pointers to the queue-complete experimental index."""
import csv,collections
with open('outputs/stage5_variant_phenotype_experimental_index.csv') as f:base=list(csv.DictReader(f));fields=list(base[0])
with open('outputs/stage5_cascade_query_coverage.csv') as f:coverage={r['variant_GRCh38']:r for r in csv.DictReader(f)}
with open('outputs/stage5_cascade_model_eqtl_comparison.csv') as f:comparison=collections.Counter(r['variant_GRCh38'] for r in csv.DictReader(f))
with open('outputs/stage5_cascade_credible_qtl_by_cell.csv') as f:credible=collections.Counter(r['variant_GRCh38'] for r in csv.DictReader(f))
assert len(base)==488 and {r['variant_GRCh38_ref_alt'] for r in base}==set(coverage)
extra=['cascade_query_status','cascade_partial','cascade_cell_contexts_returned','cascade_credible_molecular_rows','cascade_model_eqtl_comparison_rows','cascade_interpretation']
with open('outputs/stage5_variant_phenotype_experimental_multiome_index.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields+extra);w.writeheader()
 for r in base:
  k=r['variant_GRCh38_ref_alt'];c=coverage[k];r.update(cascade_query_status=c['status'],cascade_partial=c['partial'],cascade_cell_contexts_returned=c['n_cell_contexts'],cascade_credible_molecular_rows=credible[k],cascade_model_eqtl_comparison_rows=comparison[k],cascade_interpretation='Browser top-per-cell and credible subset only; no causal or independent disease claim; positional no-result is not a zero effect');w.writerow(r)
print('joined variants',len(base),'variants with model eQTL rows',len(comparison),'credible hits',len(credible))

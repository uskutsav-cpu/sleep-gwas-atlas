# GTEx v8 brain lead-eQTL lookup — exploratory

This is an exact lead-coordinate lookup for the protected 27 lead rows (26 unique variants) in two
prespecified eQTL Catalogue GTEx v8 datasets: brain cortex and brain frontal cortex (DLPFC).
The official release-7 dataset metadata, GRCh37-to-GRCh38 chain, exact rsIDs, and unordered SNP alleles
are bound in the provenance receipt. All nominal associations at matching lead coordinates are retained.

Queries: 52; lead-by-tissue rows: 54; exact matched association rows: 826.
Statuses: {'EXACT_LEAD_QTL_ROWS_PRESENT': 46, 'NO_EXACT_LEAD_QTL_ROW': 8}.

**EXPLORATORY_LEAD_LEVEL_ONLY.** These rows do not establish locus-wide eQTL support, a significant
QTL after gene/tissue multiplicity correction, colocalization, a causal gene, or tissue specificity.
A missing exact lead row does not exclude an eQTL elsewhere in the locus. Splicing QTL and full-region
molecular colocalization remain unrun; the Catalogue leafcutter folders inspected here expose selected
conditional association files but not the full nominal all-association file required for unbiased coloc.
The canonical LAVA decision and final tiers remain BLOCKED_LAVA.

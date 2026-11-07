# Brain6 evidence graph

This machine-readable graph organizes the existing selected evidence for inspection. It contains **429 variant nodes, six candidate-region nodes, 23 source-specific trait nodes and 17 gene nodes**, with 1,167 typed edges. The four extra variants come from the six post-selected frailty candidates; two of those six already belong to the Brain6 set. The graph is not an independently validated model, a causal interaction network, or a clinical tool. The original 425 variants are the selected single-signal credible-set union; the six added frailty records were post-selected. Neither is a genome-wide variant universe.

| Edge type | Meaning and limits |
|---|---|
| `in_selected_region` | Variant belongs to a selected Brain6 region. Region coordinates are GRCh37; variant node labels show the checked GRCh38 representation. |
| `conditional_abf_credible_set_member` | Variant is in the reported trait-specific set; `value` is its conditional single-signal PIP and `detail` records source effect allele. Missing edge is not proof of no association. |
| `positional_gene_overlap` | A selected representative overlaps an Ensembl gene in GRCh38. It does not establish a regulated gene, an effect, or causality. Only the subset queried in the source context table is represented. |
| `gtex_v8_significant_eqtl` | An exact-allele GTEx significant-pair lookup returned an expression association. `value` is the nominal pair p-value; details retain tissue, GTEx gene name and versioned GENCODE ID, alternate-allele NES and response hash. It is not GWAS–QTL colocalization. The graph excludes sQTL records; see the original QTL table for them. |
| `ldsc_genetic_correlation` | Genome-wide common-variant LDSC r_g and its reported SE for a specified trait pair and source. This is not a single-variant effect, local correlation, causal estimate or individual risk. Source inputs, supports and overlap differ; multiplicity corrections remain in the linked atlas rather than the graph edge. |
| `postselected_marginal_association` | Marginal p-value from the MVP insomnia, FinnGen HFRS or Ye FFS lookup of a post-selected frailty candidate. It is not a corrected p-value or independent replication. Source allele conventions and the assumed Song build require review. |
| `nearest_positional_gene` | Nearest-gene annotation and distance from the frailty source; neither establishes a causal gene or mechanism. |
| `reference_ld_r2_undirected` | An observed pairwise reference LD value from 503 1000 Genomes EUR participants, with counted alleles and signed r retained. It is not a cohort-specific causal link. |

The 100 genome-wide edges comprise 96 joint-matrix cells, UKB-insomnia–HFRS, registry-insomnia–FI and –FFS, and the native FI–FFS estimate. These are sourced from different analyses and should not be treated as one homogeneous or independent correlation matrix.

An edge is a source record rather than a validated mechanistic relationship. Nodes use rsIDs for readability and retain GRCh38 REF/ALT in labels, while eQTL joins explicitly verify the exact GTEx variant representation. Some rsIDs can have other alleles outside this graph. Ensembl and GTEx gene names can differ even when GENCODE IDs match; both source names are preserved in eQTL details. No direction from AlphaGenome ALT/REF is interpreted as a disease-risk direction.

`nodes.tsv` and `edges.tsv` are simple typed tables, and `receipt.json` records counts and source hashes. Run `python3 outputs/shared-genetics-atlas/evidence_graph/build.py` from the workspace root. The graph contains no AlphaGenome numeric score, but some gene and GTEx queries were **selected in part using AlphaGenome results**. Do not assume that this graph or any derived features are permitted for model training: check the governing AlphaGenome output terms for the intended use. This graph was not used to train any ML model. A legitimate benchmark would require independently ascertained labels, block- or chromosome-held-out evaluation, appropriate baselines and a review of terms for every input; graph paths alone do not establish mechanism or predictive performance.

[Variant trace guide](VARIANT_TRACE_GUIDE.md) demonstrates how to inspect rs62264764 and a one-LD-hop connection without treating the path as a causal claim.

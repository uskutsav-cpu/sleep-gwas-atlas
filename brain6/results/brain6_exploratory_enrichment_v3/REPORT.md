# Brain6 exploratory positional enrichment: matched 18-region sensitivity

**Result: the separate 18-region sensitivity ran successfully; no tissue or pathway passed the prespecified 5% false-discovery threshold.** This is a limited exploratory positional analysis. The original 20-region enrichment remains **NOT_ESTIMATED**, as documented in the [full-set v2 audit](../brain6_exploratory_enrichment_v2/REPORT.md). Nothing here changes a PLACO candidate, LAVA validity, QTL colocalization, replication status, or causal-gene assignment.

## Frozen scope and matched null

The [v3 protocol](PROTOCOL.md) and complete [source binding](frozen_binding.json) were committed as `ef208af9` before any enrichment statistic or p-value was calculated. The subset contains **18/20 geographic regions and 21/25 pair-specific PLACO candidates**. The other two regions (four pair-specific entries) failed the v2 pre-outcome minimum of 10 eligible matched controls: chr3:52,219,816–53,219,816 had four and chr17:43,049,526–45,343,136 had seven. Their [all-20 status](matching_status_20.tsv) remains explicit. The 18 retained regions had **861 total eligible anchor windows**, 11–103 per region (median 42.5), under the unchanged exact-width, same-chromosome, LAVA-block-span, pair-specific PLACO tested-density, gene-count, and all-20 candidate-exclusion rules. No threshold was relaxed.

The fixed 18-region subset overlaps **116 distinct GENCODE v19 protein-coding gene bodies**, 113 of which map to GTEx v8 median TPM. **17/18** regions overlap a protein-coding gene body; the gene-empty region contributes zero to the tissue statistic. These are positional overlaps, without effector-gene or causal interpretation.

With seed 20260928, **10,000** independently drawn, within-replicate nonoverlapping matched null sets were accepted in **11,082** attempts. Every accepted replicate contained one eligible null window for each of the 18 regions. The [run receipt](run_summary.json) records the gate result. Empirical one-sided high-tail p-values used `(1 + #null >= observed) / 10001`; Benjamini–Hochberg correction was separate across the frozen families of [54 GTEx bulk tissues](FROZEN_TISSUES_54.tsv) and [1,680 Reactome pathways](FROZEN_PATHWAYS_1680.tsv).

## Tissue and pathway findings

The [complete 54-row tissue table](tissue_bulk_gtex_v8_enrichment_18.tsv) has **0/54 at q < 0.05**. The smallest nominal p is 0.00210 for cervix endocervix, with q=0.11339. Among the 13 bulk brain tissues, the smallest nominal p is 0.29777 for frontal cortex BA9, with q=0.89331. These bulk-tissue results are not cell-type enrichment.

The [complete 1,680-row pathway table](reactome_positional_enrichment_18.tsv) has **0/1,680 at q < 0.05**. The smallest nominal p is 0.000400 for Golgi-to-ER retrograde transport (four region hits), with q=0.24264. The smallest adjusted q is 0.24264 across several related Reactome pathways. Reactome terms are nested and gene-body assignment is positional, so even a nominally small p is not an effector-gene or mechanism finding.

| Modality | Fixed family | Tested in v3 | q < 0.05 | Smallest p | Smallest q |
|---|---:|---:|---:|---:|---:|
| GTEx v8 bulk tissue specificity | 54 | 54 | 0 | 0.00210 | 0.11339 |
| Reactome positional pathway occupancy | 1,680 | 1,680 | 0 | 0.000400 | 0.24264 |
| Cell types | no cell-resolved panel | 0 | not applicable | — | — |
| Regulatory features | no matched genome-wide feature panel | 0 | not applicable | — | — |

The two omitted regions have unusually sparse matched controls, and their omission is related to gene density and tested-variant structure. The 18-region p-values cannot be extrapolated to all 20 regions or all 25 pair-specific candidates. Absence of adjusted enrichment in this subset is not evidence that the omitted loci or any tissue/pathway lack biological roles. No independent two-trait locus replication or qualifying LAVA result is supplied by this analysis.

## Reproduction

Run `python3 brain6/scripts/run_exploratory_enrichment_v3.py enrich` from the repository root with the mounted SSD source files identified by the [v2 receipts](../brain6_exploratory_enrichment_v2/source_receipts.tsv). The script verifies the committed v2 derived-input and v3 manifest SHA-256 values before computing outcomes. Its [output provenance](provenance.json) records the resulting table hashes. The independent v2 preparation command and complete 20-row matching audit remain [here](../brain6_exploratory_enrichment_v2/matching_feasibility.tsv).

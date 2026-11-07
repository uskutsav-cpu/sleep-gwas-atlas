# Brain6 × AlphaGenome: measured brain eQTL context

The tables join the six protected candidate IDs to the existing exploratory GTEx v8 exact-lead eQTL lookups. These are observed molecular associations from two bulk brain datasets; they are **not AlphaGenome tissue-specific predictions**. A lead is joined to AVI only when its candidate ID, rsID, GRCh38 coordinate and unordered alleles match. Leads outside the scored credible-set union are retained with an explicit missing AVI status.

## Chromosome 5: rs2431108

The existing Brain6 lookups contain two nominal associations for the same molecular trait, ENSG00000283462, at chr5:104612267:T>C (GRCh38). Both have the C allele as the QTL effect allele.

| GTEx tissue | Metadata sample size | Nominal eQTL p-value | Beta for C | Standard error |
|---|---:|---:|---:|---:|
| brain_cortex | 205 | 0.291197 | -0.0590737 | 0.0558119 |
| brain_frontal_cortex | 175 | 0.0989162 | 0.120219 | 0.0724301 |

The signs differ across these two estimates, and their uncertainty is substantial. Neither result establishes a significant molecular effect after multiple-testing correction, a tissue-specific effect, or a causal gene. A non-significant nominal result also does not exclude an effect. These are lead-only observations; they do not substitute for a full-region analysis.

Brain6’s final component-GWAS–eQTL status for this candidate is `NOT_ESTIMATED_NO_ELIGIBLE_EQTL_CONTEXT`, and pair-level GWAS–QTL colocalization is `NOT_ESTIMATED`. The exact-lead rows did not satisfy the separate eligibility gates for regional colocalization. AVI scores, a positional gene and nominal eQTL associations cannot override those gates or the failed primary LAVA family.

## Scope and files

The two original query tables have 12 lead-by-tissue status rows and 22 exact-lead gene-association rows across the six regions. `lead_eqtl_context.tsv` contains all matching gene rows, while `lead_eqtl_query_status.tsv` keeps the complete query statuses, including no-row results. Sample sizes are source metadata and are not necessarily per-variant analyzed N.

The eQTL lookup was exploratory and did not apply a gene/tissue multiple-testing correction. A missing lead-level row does not mean no regional eQTL exists. QTL beta is expressed for the listed QTL ALT allele and must not be compared with GWAS beta without checking orientation. Existing model limitations and lack of independent two-trait locus replication are unchanged.

[Source project eQTL report](https://github.com/uskutsav-cpu/sleep-gwas-atlas/blob/303c905c229104b7bf623df6de7431abb87e0727/brain6/results/exploratory_five_track_v4/lead_eqtl_context/report.md) · [Brain6 final results](https://github.com/uskutsav-cpu/sleep-gwas-atlas/blob/303c905c229104b7bf623df6de7431abb87e0727/brain6/paper/final_package_v1/BRAIN6_FINAL_RESULTS.md) · [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/)

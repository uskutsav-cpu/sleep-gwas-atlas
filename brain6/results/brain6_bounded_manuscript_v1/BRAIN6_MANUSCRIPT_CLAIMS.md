# Manuscript claim ledger

Each sentence below has a status for use in manuscript drafts. The status applies to the wording as written.

| ID | Status | Wording and evidence |
|---|---|---|
| C01 | SUPPORTED_NOW | The canonical seven-trait v3 LAVA family evaluated all 17,465 trait–locus cells, with 13,745 `TESTED`, 3,720 `NOT_RUN`, and no numerical `FAILED` cells. [Frozen status](../lava_longsleep_source_rescue_v1/BRAIN6_CONFIRMATORY_STATUS.md). |
| C02 | SUPPORTED_NOW | The family exceeded its frozen 873-cell `NOT_RUN` ceiling and was `FAILED_QC_NOT_PROMOTED`; no local-rg or final shared-locus tier was promoted. [Frozen status](../lava_longsleep_source_rescue_v1/BRAIN6_CONFIRMATORY_STATUS.md). |
| C03 | SUPPORTED_NOW | The five locked PLACO screens retain 25 pair-specific candidate loci grouped into 20 geographic regions. [Candidate table](candidate_evidence_25.tsv); [region table](region_evidence_20.tsv). |
| C04 | SUPPORTED_NOW | Long sleep contributed 1,291 canonical `NOT_RUN` cells, the most among the seven traits. [Trait receipt](../lava/canonical_v3_not_run_by_trait_v1.tsv). |
| C05 | SUPPORTED_NOW | The source-verified 88-locus total-N linear pilot gained zero net `TESTED` loci and failed its predeclared advancement rule. [Pilot adjudication](../lava_confirmatory_pilots_v1/longsleep_linear_88_terminal_20260927.md). |
| C06 | EXPLORATORY_ONLY | FinnGen provides global pair-level directional replication for insomnia–ADHD; it does not independently replicate any candidate locus. [Replication table](replication_25.tsv). |
| C07 | EXPLORATORY_ONLY | A phenotype-adjacent Yale ≥10-hour UKB+MVP file gave three exact-variant P/direction labels among the 25 candidates; its UK Biobank overlap prevents independent replication. [Replication table](replication_25.tsv). |
| C08 | EXPLORATORY_ONLY | Candidate variants in six pair rows intersect published GTEx brain molecular credible sets; their reported PIP values are molecular QTL PIP, not GWAS PIP. [QTL table](qtl_context_25.tsv). |
| C09 | EXPLORATORY_ONLY | Six Ensembl regulatory features overlap five exact leads; overlap is coordinate context and does not assign causal genes. [Regulatory table](regulatory_gene_25.tsv). |
| C10 | REQUIRES_LAVA_CONFIRMATION | Any candidate demonstrates a shared local genetic correlation or qualifies for a frozen shared-locus evidence tier. [Frozen tier rule](../../config/shared_locus_evidence_tiers_v1.json). |
| C11 | REQUIRES_LAVA_CONFIRMATION | A candidate has GWAS fine-mapped causal-variant PIP, trait–trait or GWAS–QTL colocalization, or enriched tissue/cell type/pathway. These claims also need their separate method inputs and validation; [method table](method_feasibility_25.tsv). |
| C12 | REQUIRES_LAVA_CONFIRMATION | Any of the 25 candidate loci replicates as a two-trait shared signal in a nonoverlapping cohort. This additionally needs independent two-trait data; [replication table](replication_25.tsv). |

Do not replace `NOT_ESTIMATED` with a negative biological conclusion. Any future change in claim status requires a new source-bound analysis and should not edit these frozen descriptive receipts in place.

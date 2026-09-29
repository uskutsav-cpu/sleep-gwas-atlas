# Brain6 final package readiness

**Disposition: complete transparent secondary-validation paper (Success B).** The original primary LAVA endpoint remains `FAILED_QC_NOT_PROMOTED`; this package does not represent a LAVA rescue or a confirmed shared-locus manuscript. It reports all prospectively executable secondary and exploratory analyses with explicit source and method holds. No protected candidate has a positive A–D evidence class.

## Delivered scope

| Component | Final state |
|---|---|
| Global atlas and PLACO | Locked 72-row global display; 25 pair-specific PLACO candidates in 20 geographic regions; original correction and identities preserved. |
| Primary local test | 17,465 canonical LAVA cells; 13,745 `TESTED`, 3,720 `NOT_RUN`, 0 `FAILED`; frozen ceiling 873; no promotion. |
| Source rescue | Exact long-sleep public/local search and LAVA method-contract audit complete; all six other trait source ledgers complete. Insomnia zero-row engineering fix and 88-locus pilot complete, but its advance gate failed. Validated family recovery zero. |
| Secondary local method | Prospectively selected SUPERGNOVA; two eligible pairs complete with four workers each and all 1,693 blocks per pair; three long-sleep pairs prospectively inapplicable; 8,465 family slots and 25/20 crosswalks. Zero protected candidates supported. |
| GWAS fine mapping and trait coloc | SuSiE LD gate failed explicitly; separately frozen ABF single-signal method completed 28/50 trait units and six pair coloc tests, with all holds and prior sensitivity. |
| QTL and replication | Four actual ADHD-component brain GWAS–QTL tests; no pair-level QTL support. All 25 candidate replication statuses explicit; none independent two-trait locus replication. |
| Regulatory and enrichment | Positional molecular, Ensembl, gene, and literature context retained. Full 20-region enrichment failed matched-control gate; frozen 18-region sensitivity tested all 54 tissues and 1,680 pathways with no FDR discovery. No reviewed cell-resolved or regulatory-feature enrichment universe. |
| Integration | Seventeen SHA-256-bound tables, six manuscript documents, 25 candidate rows, 20 geographic rows, full primary and secondary family tables, and independent final audit. |

## Reproduction and decision record

The exact [source-to-output hash manifest](BRAIN6_FINAL_PROVENANCE.json) records frozen anchors, all source files, each of the 17 tabular outputs, and the builder hash. The [independent audit receipt](BRAIN6_FINAL_AUDIT.json) records document hashes and final cross-table checks. Source URLs and release identifiers are retained in the [GWAS availability manifest](../../manifests/gwas_external_availability.tsv), [QTL input receipts](../../results/brain6_exploratory_functional_v1/coloc_input_receipts.json), [enrichment source receipts](../../results/brain6_exploratory_enrichment_v2/source_receipts.tsv), and method reports. Software versions and exact scripts/patches are bound by their respective upstream receipts; `coloc` was 5.2.3 and the SUPERGNOVA implementation/patch is hashed in its protocol.

From the repository root, `python3 -m unittest brain6.scripts.tests.test_build_brain6_final_package_v1 -v` exercises all 25/20 joins and final table shapes using synthetic secondary statuses; `python3 brain6/scripts/audit_brain6_final_package_v1.py` independently checks the real package. The builder refuses an existing final directory, and the audit receipt is write-once. Current output directories and prior LAVA v2/roundoff/pilot receipts remain preserved.

The strongest manuscript wording is in [results](BRAIN6_FINAL_RESULTS.md) and the [claim ledger](BRAIN6_FINAL_CLAIMS.md). The [discussion](BRAIN6_FINAL_DISCUSSION.md) and [limitations](BRAIN6_FINAL_LIMITATIONS.md) distinguish numerical non-support, method inapplicability, and genuine missing source information. A future confirmatory analysis requires new exact-source information and a new versioned, prospectively frozen complete-family workflow; no current table should be edited in place to imply promotion.

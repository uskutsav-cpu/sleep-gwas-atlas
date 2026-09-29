# Brain6 bounded manuscript readiness and handoff

## Finished and verifiable now

| Work | Status | Source/output |
|---|---|---|
| Canonical seven-trait LAVA v3 univariate family | COMPLETE, FAILED_QC_NOT_PROMOTED | 17,465/17,465 cells; 13,745 `TESTED`, 3,720 `NOT_RUN`, 0 `FAILED`; frozen ceiling 873. [Terminal status](../lava_longsleep_source_rescue_v1/BRAIN6_CONFIRMATORY_STATUS.md). |
| Pair-specific candidate identities and geography | COMPLETE EXPLORATORY | 25 pair candidates, 20 geographic groups; [candidate table](candidate_evidence_25.tsv), [region table](region_evidence_20.tsv). |
| Independent genotype LD and signed-direction diagnostics | COMPLETE WITH EXCEPTIONS | Five cross-pair edges; six out-of-range candidate raw signed correlations and 24 missing exact reference alleles retained in [v6](../exploratory_five_track_v6/region_evidence_report.md). |
| Global pair and phenotype-adjacent replication classification | COMPLETE EXPLORATORY | Six insomnia–ADHD rows carry **global** FinnGen context; Yale long-sleep candidate lookups classified 3/14/2/6; [25-row table](replication_25.tsv). |
| Published molecular, positional, regulatory, prior-literature context | COMPLETE DESCRIPTIVE | [QTL](qtl_context_25.tsv), [tissue](tissue_context_25.tsv), [regulatory/gene](regulatory_gene_25.tsv), [24 prior-study records](prior_literature_records.tsv). |
| Explicit fine-map, coloc, and pathway status tables | COMPLETE READINESS CLASSIFICATION | [GWAS fine-mapping](fine_mapping_25.tsv), [trait–trait coloc](trait_trait_coloc_25.tsv), [eQTL/sQTL](eqtl_sqtl_evidence_25.tsv), [pathway](pathway_20.tsv); all missing method prerequisites remain visible. |
| Manuscript drafting | READY FOR BOUNDED EXPLORATORY REPORT | [Results](BRAIN6_RESULTS.md), [Methods](BRAIN6_METHODS.md), [Discussion](BRAIN6_DISCUSSION.md), [Limitations](BRAIN6_LIMITATIONS.md), [claim ledger](BRAIN6_MANUSCRIPT_CLAIMS.md). |

## Analyses not estimable from currently verified inputs

| Analysis | Specific missing prerequisite | State |
|---|---|---|
| Independent locus-specific two-trait replication | Nonoverlapping cohort with aligned GWAS for both traits at candidate loci | NOT_ESTABLISHED; Yale includes UKB and differs in long-sleep cutoff; SleepChart variant data inaccessible without approved access. |
| GWAS fine-mapping | Reviewed trait-specific N/binary semantics, full allele-aligned regional signed LD, pinned SuSiE runtime | NOT_ESTIMABLE; no GWAS PIP. |
| Trait–trait colocalization | Reviewed paired regional data, N, overlap/prior contract, suitable runtime; multi-signal path requires fine-mapping | NOT_ESTIMABLE; no PP.H4. |
| GWAS–eQTL/sQTL colocalization | Source-bound full-cis QTL summary data with harmonized alleles/build/N and priors | NOT_ESTIMABLE; molecular credible sets are descriptive. |
| Tissue/cell-type or pathway enrichment | Reviewed locus-to-gene weights, tissue/cell background, gene universe, pathway set and multiplicity family | NOT_ESTIMABLE; source tissue labels only. |
| Full bivariate LAVA, final protected tiers | New source-valid seven-trait family passes all frozen gates and reviewed pair overlap | BLOCKED_LAVA; no candidate promoted. |

## Exact external information needed

The smallest path is the response described in the [unsent custodian request](../lava_longsleep_source_rescue_v1/LONG_SLEEP_DATA_REQUEST.md): a variant-keyed analyzed `N`/`NMISS` (ideally per-variant cases and controls) for the exact Dashti ≥9 h versus 7–8 h release, plus a written statement binding the released BETA/SE/P to the original association model and any transformation. The source must define effect allele/build, missingness, case fraction or continuous-model mapping, and access rights. A new exact-phenotype rerun may substitute only if its phenotype and full per-variant statistical semantics are source-verified. A ≥10 h or continuous-duration phenotype remains a separate exploratory/sensitivity phenotype. Other trait repairs and pair-overlap re-estimation remain necessary because perfect long-sleep repair alone leaves 2,429 `NOT_RUN` cells.

## Restart and readiness decision

The exact initial command is in [CONFIRMATORY_RESTART.md](CONFIRMATORY_RESTART.md). The manifest template currently returns `BLOCKED`; when source and signed method review arrive, the guard validates frozen anchors and new hashes, then freezes an isolated admission receipt. A prospective 88-locus pilot, multi-trait projection, new seven-trait family, frozen QC, bivariate gate, and 25-candidate rejoin follow in that order. The old v3 runner cannot be repointed to a new source because its identity is hardcoded.

**Manuscript readiness:** the descriptive atlas, transparent canonical failure, and limitations are ready for a bounded exploratory manuscript or supplement. The project is **not ready for a confirmatory shared-locus or causal-mechanism manuscript claim**. This readiness decision is tied to method inputs and frozen QC rather than the direction or significance of candidate associations.

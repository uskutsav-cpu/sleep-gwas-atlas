# Source eligibility, access and actual phenotype sensitivity

The public FinnGen R13 browser supports a bounded clinical phenotype sensitivity at three fixed biallelic variants. It does not provide an eligible complete dense regional source in this audit. No fully independent, source-eligible two-trait replication pair has been acquired.

The official [access page](https://www.finngen.fi/en/access_results) explicitly distinguishes public browsing from full summary-statistic cloud access requiring a form. Three ordinary browser API objects were acquired after a frozen manifest, each <1MB. No form, restricted cloud file, authentication or access override was used. Bulk downloads and individual participant data were not attempted. Cite FinnGen/Kurki2023 and acknowledge participants/investigators in any resulting publication. A public webpage's accessibility does not imply arbitrary raw GWAS redistribution rights.

R13 was publicly released 2June2026, with 500,186 individuals and 21,311,644 variants. The official release has 2,755 endpoints; the current `/api/phenos` response has 2,754 rows and the deployed browser introductory text mentions an older 2,466 binary/3 quantitative count. These discrepancies are preserved. Clinical counts match between current phenotype metadata and selected native rows: insomnia51,643/446,273; ADHD5,559/489,493; depression61,203/432,672 cases/controls. Counts are endpoint metadata, not verified variant-level observed sample sizes; the annotation NS519,972 is not substituted.

All browser positions are GRCh38. Official effect/AF convention is ALT; `maf` is ALT frequency despite its name. R13 uses REGENIE v3.3 with approximate Firth/LOCO, age, sex, ten PCs, chip version and legacy batch. The frozen protocol's generic SAIGE reference is explicitly corrected in [processing correction](FIXED_VARIANT_PROCESSING_CORRECTION.md); the immutable original and all frozen statistical rules remain intact. Native P is authoritative; reconstructed Wald P is diagnostic.

The advertised `/api/variant/chr-pos-ref-alt` route returns the exact variant identity, annotations and phenotype rows. The public backend inserts null-statistic placeholders for phenotype rows filtered from its display matrix. Actual maximum nonnull P across these objects is near .05, consistent with display censoring; the deployed threshold is unverified. Therefore an absent/native-null association is non-estimable rather than a null or P=1. The region metadata route does not return dense associations; the association route also applies a deployed unknown `p_threshold`. No dense regional colocalization is admissible from this acquisition.

The frozen four-slot intersection–union test uses `max(native P_sleep,native P_outcome)` and .05/4=.0125. All four hypotheses stay in the denominator. Corrected outputs give:

| Fixed hypothesis | Conjunction P | Bonferroni P | Decision |
|---|---:|---:|---|
| A rs7105462 insomnia–ADHD | unavailable | unavailable | NOT_ESTIMATED: ADHD placeholder |
| A rs7105462 insomnia–MDD | .00275816346 | .01103265382 | Both clinical associations at fixed variant |
| B rs9485410 insomnia–MDD | .03164171394 | .12656685578 | Did not meet family threshold |
| C rs77960 insomnia–ADHD | .000346760803 | .00138704321 | Both clinical associations at fixed variant |

Six eligible rows have finite beta/SE/native P/ALT AF; three exact variant identities match the frozen manifest. The parent independently validates the raw JSON, identities and calculations in [independent checks](independent_raw_JSON_checks.tsv). The initial parser's conservative holds are preserved separately; the correction only recognizes native `chr` rather than expected `chrom`. The original manifest/hash and hypotheses were not changed.

At A ALT A has negative effects for both insomnia and depression; at C ALT A has positive effects for both insomnia and ADHD. These are within-FinnGen descriptive directions, not a frozen discovery-concordance endpoint. The intersection–union test supports both clinical associations at a fixed variant. It does not establish that their causal variants are shared, nor regional covariance, new locus discovery, disease mechanism or fully independent replication.

## Cohort graph and resolved ADHD roster

[Source graph](SOURCE_COHORT_GRAPH.tsv) gives 32 typed edges, separating contributing cohorts, biobank organisations, endpoint subsets and legacy annotation labels. New official Demontis2023 Supplementary Table1 and Figshare ADHD2022.xlsx independently agree on thirteen cohorts and all26 case/control counts, total38,691/186,843. The ten PGC cohorts are Barcelona, Bergen, Cardiff, CHOP, Germany, IMAGE-I, IMAGE-II, PUWMa, Toronto and Yale-Penn; the other components are iPSYCH1/2 and deCODE. The roster includes age, sex fractions, covariates and case-control versus trio designs; trio controls must not be treated as ordinary unrelated population controls.

No exact PGC roster name matches the ten declared FinnGen biobank organisation names. These naming units differ, and participant linkage, legacy control-pool aliases and endpoint-specific cohort allocation are unavailable. Public FinnGen point annotations carry NFBC66/NFBC86 legacy AF/INFO dataset labels; such annotation is not evidence of those cohorts' endpoint case/control membership. Absence of NFBC from the thirteen ADHD names is informative source documentation, not proof that no individual/control is reused. [Overlap decisions](COHORT_OVERLAP_DECISIONS.tsv) retain UNVERIFIED for discovery–FinnGen. The new four-slot experiment is PHENOTYPE_SENSITIVITY because registry diagnoses differ from questionnaire sleep complaints and meta-analysis case definitions.

The official Demontis [author correction](https://doi.org/10.1038/s41588-023-01350-w),1March2023, fixes an author first name and adds a 23andMe acknowledgement. It does not report a change to these thirteen GWAS cohort counts. An acknowledgement alone is not added as a contributing GWAS cohort. The correction source is acquired and hash-recorded with the other metadata.

## Independent pair feasibility

Partners Biobank primary data availability explicitly restricts underlying data for IRB reasons and offers summary statistics on corresponding-author request. Its European patient cohort has3,135 cases/14,920 controls and5,508,534 QC-passed autosomal biallelic SNPs, but no authorized summary object was acquired. STARRS gives a distinct lifetime-insomnia cohort:3,237/14,414 all-ancestry cases/controls, EUR cohort totals4,756/1,817/4,900, explicitly GRCh37/hg19. Exact EUR case/control counts must not be inferred by rounding the published prevalence percentages; no complete authorized compatible summary object was acquired. Source coverage and discovery/control overlap remain unresolved. HUNT-only source release and rights remain unverified. Pairing any of these with FinnGen ADHD/MDD is a potential route, not an admitted independent pair. [REPLICATION_MASTER.tsv](../REPLICATION_MASTER.tsv) contains all actual and prospective options; no row is promoted to INDEPENDENT_TWO_TRAIT.

Original insomnia reused with a different psychiatric GWAS is one-trait evidence at most; independence of that psychiatric source still requires cohort auditing. PGC2025 noUKB also reuses older PGC cohort names in the historical audit, so changing release is not itself independence. The new ADHD source metadata/readme explicitly prohibit publicly reposting GWAS result files despite generic Figshare CC BY metadata. Raw discovery subsets and raw source objects remain ignored; source facts and derived analytical results carry exact hashes and terms.

[Source feasibility](SOURCE_FEASIBILITY.tsv), [public acquisition ledger](PUBLIC_ACQUISITION_LEDGER.tsv), [investigator request drafts](INVESTIGATOR_REQUESTS_NOT_SENT.md), and [unresolved literature access](../UNRESOLVED_LITERATURE_ACCESS.tsv) identify the executable next actions and genuine external requirements. No requests were sent.

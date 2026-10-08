# Sleep-specific novelty and prior-art audit

Execution/access date: **2026-10-08**. This is technical evidence documentation, not manuscript text. The audit is bounded, and neither an exhaustive systematic review nor a priority certification. Search indexing may lag publication; the execution date is not proof that every publication through that date was indexed.

**Decision: no adequately supported novel-priority result is established.** Direct review of prior result tables materially weakens the historical novelty labels. The reproducible receipt is `logs/literature_audit_receipt.json`; every original pair is retained in `tables/NOVELTY_MASTER_1200.tsv`, all historically replicated pairs in `tables/NOVELTY_REPLICATED_23.tsv`, and concrete quantitative comparators in `tables/novelty_crosswalk.tsv`.

## Executed evidence work

- Re-extracted and compared all **8,418** Morrison 2024 LDSC rows against the original retained XLSX: phenotype labels and rg/SE/P values agree with the six historical extracted sheets, within floating-point representation tolerance.
- Inspected Goodman 2025 SD24 directly: **380** nonblank phenotype rows. The paper describes 375 representative traits and 2,250 tests. The 380 rows must not be reported as 380 independent external discoveries or assumed to mean a 2,280-test declared family.
- Inspected Dashti 2019 Data18: **45 displayed trait rows**; absence from this selected table does not establish absence from the broader 224-trait screening family. Data17 is a TWAS table, and is explicitly excluded as rg prior evidence.
- Inspected Li/Zhao 2020 original DOCX table: **53 data rows** after two header rows, with device-derived constructs. Different device sleep cutoffs must not be equated to the self-report short/long sleep definitions.
- Retrieved and inspected the complete Nature 2026 Sleep Chart workbook, SF5: **1,054** original-study comparisons; **309** reported Dashti-comparator rows; **361** Austin-Zimmerman-comparator rows. The latter two are selected displayed results and cannot certify that undisplayed pairs were never studied. The [30July2026 author correction](https://www.nature.com/articles/s41586-026-10920-x) adds two omitted prior references; it does not state a numerical SF5 correction.
- Retrieved current Communications Medicine 2026 disease correlation Supplementary Data11: **4,046** displayed rows, together with sleep/data-source metadata and related supplementary files. The publication reports 34 sleep GWAS and 113 disease datasets; additional BMI-adjusted constructs and incomplete rows explain why a simple 34×113 multiplication does not equal the displayed total. This audit does not alter the prior paper's family.
- Retrieved and inspected original Wang2019 sleepiness Data6 (232 named rows; paper text states233 traits while its figure refers to224 available) and Dashti2021 napping Data10 (257 rows; both BMI-adjusted/unadjusted models). Birth-weight comparators are recorded; smoking labels are not equated to tobacco disorder.
- Retrieved Sun 2022 full primary PDF and programmatically extracted all ten relevant sleep–GERD estimates from Table2, including the FinnGen R6 outcome-side comparisons.
- Ran focused live primary-source literature searches for broad sleep phenome screens, psychiatric genetics, GI, pain, smoking, COPD, musculoskeletal outcomes, and unresolved replicated relationships. Raw query results/failed fetches are retained under `sources/literature*`.

The historical audit searched/classified 603 significant pairs on 2026-08-29; the remaining **597** pairs were not historically pair-searched. All 1,200 now have local supplement assessment and explicit search limitations. Current pair-specific database searches have **not** been rerun exhaustively for all 1,200. `NOT_SEARCHED` and `UNRESOLVED` are preserved rather than reclassified as negative evidence.

## Material correction to historical interpretation

The prior registry calls all six Morrison constructs latent sleep factors. The primary paper distinguishes three multivariate factors from **three single-indicator traits** (non-insomnia, duration, regularity). The S14 insomnia and S13 duration sheets are therefore substantially closer comparators than that registry suggests. Insomnia was reverse-coded toward better sleep health; negative prior rg values must be reoriented before comparing with this study's positive insomnia liability.

Morrison already reports insomnia with ICD K21 reflux, K44 diaphragmatic hernia, J44 COPD, and M51 disk disorders. Their retained original rg values are respectively −0.4383, −0.4641, −0.2763, and −0.3064. This reverses their sign for disease-liability interpretation, not their SE or P. These findings are **previously published substantially similar relationships**, with exact source subsets, cohort sizes, disease definitions and release differences still qualified. The old `NO_DIRECT_RG_FOUND` label for insomnia–reflux is not defensible as a current priority conclusion. Historical files remain unchanged.

Sleep Chart SF5 block b explicitly uses **Dashti sleep GWAS**, and includes short/long sleep with FinnGen reflux, COPD and intervertebral disorders. Thus a substantially similar outcome-side design already appeared before this execution. It is inappropriate to advertise these findings as first reports based only on use of FinnGen R13 instead of an earlier release. See the comparator values and locations in the crosswalk.

Sun 2022 already estimates insomnia, short/long sleep, sleepiness and napping against UKB GERD and FinnGen R6 GERD. Napping–GERD is therefore not supported as a new genetic correlation. Snoring–GERD has prior MR evidence; that is related genetic evidence, not an interchangeable rg estimate.

## Classification policy

| Requested distinction | Audit treatment |
|---|---|
| Previously published identical correlation | Requires exact sleep definition/source subset, external definition, analysis and orientation. None certified as byte-identical by this audit. |
| Published substantially similar association | Direct phenotype-level rg comparators are recorded, including differences in cohorts/releases. |
| Known clinical association without direct rg | Kept as contextual evidence only; never converted into direct rg. |
| New cohort replication of prior finding | Recorded separately from novelty; retained sleep GWAS prevents two-trait independence. Earlier FinnGen release reuse may limit even cohort-newness. |
| Potentially underreported relationship | Unresolved unless a documented breadth/synonym search plus full supplements supports it. |
| New quantitative characterization | Appropriate for updated estimates, complete null reporting and cross-source comparison; does not establish first discovery. |
| Adequately supported novel result | **Zero**. |
| Literature comparison unresolved | Explicitly retained, including source access or exact-definition uncertainty. |

The current receipt classifies **83/1,200** as having direct substantially similar published table comparators, **586** with broad external coverage under different or unresolved sleep constructs, **10** with related genetic/MR evidence, and **521** unresolved. These are conservative operational classes, not a census of all prior literature. A direct published comparator means a reported rg estimate, **not automatically a statistically established association**: prior nominal-P flags are supplied, but original family-corrected support is not recertified. Every direct comparator among the23 has at least one nominal P<0.05; duration–disk in Morrison hasP=.00583, which is not thereby a prior family-corrected result. Current significance or a different cohort cannot be turned into first-priority evidence without broader comparison. Among the 23 historical replications: **15 direct prior comparators, 6 related genetic/MR contexts, 2 unresolved**. The unresolved findings include snoring–cholelithiasis and sleep-apnea–abdominal pain; unresolved status does not support novelty.

## Reproduction and limits

Run `python discovery_extension/sleep_submission_evidence_v1/scripts/audit_literature.py` with openpyxl available. The script reads historical frozen tables and original literature workbooks, asserts the 1,200/217/603/23 universe boundaries, verifies original supplement extractions, and regenerates the tables. It does not replace historical outputs or silently launch searches. Source hashes are in `sources/literature_source_hashes.tsv`; download outcomes and license qualifications are in `sources/literature_access_ledger.tsv`. Local full-text/supplement cache is excluded from public release; see `sources/literature_CACHE_RELEASE_BOUNDARY.md`. Source-free verification uses `--verify-canonical`; official cache recovery/hash checking uses `--acquire`.

Remaining priority work: retrieve complete García-Marín 2021 SLEEP supplements and exact prior source metadata; resolve all coding/source-subset ambiguities; extend dated searches with synonyms and full-result reviews; and obtain sleep/statistical-genetics specialist review. An empty search result is never sufficient evidence of priority.

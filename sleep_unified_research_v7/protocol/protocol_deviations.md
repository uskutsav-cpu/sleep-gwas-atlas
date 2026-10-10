# V7 protocol deviations and access decisions

**Recorded:** 2026-10-10, before Fan SD2/SD4/SD11 row matching and before benchmark-task sampling or answer generation.

The initial Human GWAS ATLAS, Morrison, and Goodman match tables had already been generated and inspected before this record was written. Those initial counts are exploratory and will be replaced by one final build after the decisions below. No V7 solver answers have been generated.

## Independent-review correction to Human GWAS ATLAS source classes

Before the final build, independent statistical and scientific reviewers found that the initial 30 HGA class-A rows were assigned from publication PMID, approximate sample size, ancestry text, and phenotype labels. Those fields do not establish the preregistered source release, ascertainment/coding, cohort frame, effect direction, estimator/reference, or exact estimand. The independent statistical review checked the released pair table and found opposite-signed atlas/HGA estimates for 16 of the 30 candidates. The builder now classifies all 30 as **E, insufficient evidence**, pending those missing comparability details. Their source references remain available in full. This is a classification correction, not new biological evidence; V7 exports no HGA numeric result values.

All 30 HGA candidates also have Morrison-related B coverage, so union precedence leaves them B in the four-snapshot union. The corrected union is therefore A=0, B=1,183, C=0, D=366, E=88. The 361 qualified match-absent records are unchanged because none of the 30 metadata candidates met the all-snapshot qualification rule.

The final output preserves every generated candidate row reference rather than truncating long lists. Goodman SD24 contains six composite scores (SHS-ADD and SHS-PC1–PC5), not five. These corrections affect comparator coverage reporting only and do not change any V4 result.

## Fan et al. 2026 supplements included in row-level matching

The official publisher supplements SD2, SD4, and SD11 are already present on the attached Extreme SSD. Their SHA-256 values match the previously recorded acquisition hashes in the V4 source ledger, so they are reused in place; no verified file was redownloaded. SD11 contains 4,046 result rows, with source definitions in SD2 and SD4. V7 records only sheet row references, source matches, and derived classifications; it does not copy Fan rg, SE, P, FDR values or redistribute workbook contents. The publisher materials' reuse conditions still apply.

The article states a CC BY-NC-ND 4.0 license. V7's derived annotations and row references are not assumed cleared for public redistribution; human investigators must resolve the license and attribution before public deposit.

## SleepChart 2026 carried forward as audited context, not a V7 row-level comparator

The V4 benchmark already reviewed the relevant SleepChart supplement and its row-level schema; its compact audit is in `sleep_unified_research_v4/cell_genomics_readiness_v1/benchmark_matrix.tsv` and the associated report. The inspected workbook is not in the active checkout, and the instruction not to redownload a verified source is followed. SleepChart is therefore reported as a prior audited resource with partial, related duration–outcome coverage; V7 does not count it in the row-level union or classify it as wholly incommensurate. This is not a negative search and does not imply that the resource lacks comparable pairs.

## Scope and interpretation

No GWAS/LDSC calculations, historical estimates, source releases, or frozen V4/V5/V6 files are changed. Exact-source matching remains a resource-coverage assessment. A release-specific absence is not evidence of novelty, and no historical association is described as a new biological result.

## Experiment B task benchmark deferred before task sampling

The frozen protocol proposed a 60-task paired baseline-versus-atlas automated-answer benchmark. It is not executed in this completion. The V4 package already has a completed independent consumer review of the actual machine-readable candidate, including result lookup, source joins, validation qualifications, schema, and an isolated byte-for-byte rebuild. Repeating that review as another automated agent exercise would not establish human researcher usability, improve statistical validity, or resolve the current scientific novelty limit. A new AI-only task experiment could address comparative answer performance, but that evidence is not necessary to write a responsibly bounded resource paper and would not support a biological discovery claim. This is a scope deviation from the frozen V7 plan: no task frame, task answers, blinded adjudications, or task-performance figure will be reported, and the V7 paired-utility hypothesis remains untested. The V4 consumer-review evidence must not be described as a controlled comparison against baseline resources.

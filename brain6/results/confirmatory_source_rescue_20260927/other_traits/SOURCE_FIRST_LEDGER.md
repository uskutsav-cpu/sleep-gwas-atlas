# Non-long-sleep LAVA source-first rescue ledger (2026-09-27)

**Scope and status.** This is a read-only, outcome-blinded audit of the six other
canonical traits. No source has been admitted to a confirmatory replacement, no
pilot was run, and canonical v3, v2, roundoff, receipt, LD, and candidate files
were not changed. The frozen family has 3,720 `NOT_RUN` cells; 2,847 recoveries
are required to meet the 873-cell maximum. The source-independent lower bound
requires at least four traits to be repaired perfectly; the unique four-trait
set is long sleep, insomnia, Parkinson disease, and MDD. Those are arithmetic
constraints, not predicted recoveries. See
`brain6/results/lava_confirmatory_rescue_plan_v1/rescue_planning_report.md`.

## Direct-source evidence

| Priority | Trait; canonical `NOT_RUN` | Original release and relevant source fields | Current N mapping and prospective status |
|---|---|---|---|
| 2 | Insomnia; 671 | Jansen 2019 UK Biobank-only, GRCh37, binary field 1200, 109,402 cases/277,131 controls. Raw `SNP UNIQUE_ID CHR BP A1 A2 MAF OR SE P N INFO`. Official README explicitly defines `A1` effect/minor allele, `OR` odds ratio for A1, and `N` per-SNP sample size. Exact raw SHA-256 `32848cd92a6324c9048d…6a9065b`. | Harmonizer ignored literal per-SNP `N` because its binary branch preferred case/control metadata; canonical `N=313750.03595553286` is a **constant conventional effective N**, though source card says `total`. Best same-release candidate: source `N_i` plus source-cohort case fraction. **ADMISSIBLE_SENSITIVITY_ONLY for a predeclared technical pilot; not yet confirmatory-admitted.** |
| 3 | Parkinson; 636 | Nalls 2019 public restricted-cohort-excluded GRCh37 fixed-effect meta-analysis includes 18,618 UKB proxy cases. Raw header has `chromosome base_pair_location effect_allele other_allele beta standard_error effect_allele_frequency p_value N_cases N_controls` but no rsID or INFO. Exact decompressed source SHA-256 `409b9151871aafb5f6fe9acf0287fe9f5df448eda255d06b9d55b81ac8e68eeb`. | Harmonizer maps `N_cases`,`N_controls` to effective N, then LAVA combines that N with one aggregate observed fraction. The released counts vary sharply by SNP and the meta-analysis blends PD and familial proxy status. Total-count substitution requires a source-supported joint model, an allele-aware position-to-reference mapping, and new overlap covariance. GP2 2025 is hg38 and has no per-variant N in its public README. **METHOD_HOLD.** |
| 4 | MDD; 589 | Howard 2019 UKB/PGC public release without 23andMe, raw `MarkerName A1 A2 Freq LogOR StdErrLogOR P`; no SNP-level N, coordinates, or INFO. Raw SHA-256 `082a1602e405fa5e8b40917521f063444131db68e2f435e7aacfd80084852110`. | Harmonizer used fixed N_eff from cohort counts. PGC MDD2025 has per-variant `NEFF`, `NCAS`, `NCON` and richer coverage, but its cohort-wise effective N differs from effective N computed from pooled case/control counts. Its 111-locus two-arm trait-only pilot was numerically 108 `TESTED` in each arm and remains a sensitivity with model/overlap hold. **ORIGINAL_SOURCE_N_MISSING; MDD2025_SCIENTIFIC_HOLD.** |
| 5 | ADHD; 223 | Demontis 2023 European iPSYCH/deCODE/PGC meta-analysis, hg19, effect allele `A1`; raw `CHR SNP BP A1 A2 FRQ_A_38691 FRQ_U_186843 INFO OR SE P Direction Nca Nco`. Raw SHA-256 `c58a96031ec44b1edba81c91d603037a2efd03fcd946531424866e3f5f1d40c5`. | Harmonizer computed per-SNP N_eff from `Nca`/`Nco` and LAVA applied one aggregate fraction. Counts exist, but participating cohorts and case fraction vary. Requires source-based interpretation and fixed-fraction assessment before a pilot. **METHOD_HOLD, LOWER FAMILY PRIORITY.** |
| 6 | Bipolar; 168 | Mullins 2021 PGC3 European mixed-cohort meta-analysis, GRCh37, `A1` effect allele. Raw PGC metadata calls `NEFFDIV2` half effective N and describes `NCAS`/`NCON` as *effective sample size* cases/controls; also gives `IMPINFO`, `FCAS`, `FCON`. Raw SHA-256 `3bf42353fd8fcd832986c3ad122747fba458e9712556659fe73b16a91b3081dd`. | Harmonizer used `2*NEFFDIV2`, appropriately effective for LDSC, then LAVA used an observed aggregate fraction. `NCAS+NCON` is **not established as literal analyzed total**. Requires cohort-wise meta-analysis/LAVA model review. **METHOD_HOLD, LOWER FAMILY PRIORITY.** |
| 7 | Schizophrenia; 142 | Trubetskoy 2022 PGC3 European autosomal stratum, GRCh37, `A1` effect allele, `NCAS NCON NEFF`, `IMPINFO`, `FCAS`, `FCON`. Source `NEFF` is per-arm convention; the registered LDSC schema explicitly doubles it. Raw SHA-256 `dbba3a85575c99fd1c2e3497d0c7a44539ccfdbf2e69742f8bdbda879230bcd7`. | Harmonizer used `2*NEFF`, then LAVA applied aggregate fraction. `NCAS+NCON` cannot be substituted merely because the counts are present; PGC case-control/trio conventions and meta-analysis need a supported LAVA mapping. **METHOD_HOLD, LOWER FAMILY PRIORITY.** |

Paths and checksums come from `brain6/manifests/gwas_external_availability.tsv`,
`brain6/qc/source_archive_audit.tsv`, and direct inspection of each mounted
source's literal header. The insomnia gzip was independently SHA-256 checked
again on this date. The archived harmonized hashes and source cards are in
`brain6/manifests/locked_dense_input_audit.tsv` and
`extensions/brain6/work/overnight-v03/source-cards/`.

### Insomnia full-file source-only scan

The exact raw Jansen gzip yielded 10,862,567 rows; all had parseable `N`.
Minimum `N=366461`, maximum `N=386533`, mean `N=384891.81808001734`;
none was below the canonical constant `313750.03595553286`. This scan read
only the `N` field, not effects, SEs, P values, LAVA outcomes, or candidates.
For the first exact source SNP (`rs12184267`), raw `N=381202`; the canonical
harmonized row has `N=313750.03595553286`. The direct source `N` range is
close to the article's 386,533 UKB participants and varies with variant
availability. That supports—but does not independently prove—the
interpretation as per-variant analyzed participants. Per-variant case/control
counts are absent, so the case fraction is not directly observed for each SNP.

The [official CNCR release page](https://cncr.nl/research/summary_statistics/)
identifies the Jansen wave-2 public statistics and explicitly says 23andMe is
excluded. Its adjacent [exact-file README](Insomnia_sumstats_Jansenetal.readme_220525.txt),
downloaded from the same public share as the data file, is 1,300 bytes with
SHA-256 `a5e92e264da88054173fc441e05f4587a6c6134078cc12b69e7294986fcc84a0`.
It explicitly defines `N` as **per-SNP sample size**, `OR` as odds ratio for
the A1 effect allele, GRCh37 position, MAF, and INFO. This is direct source
evidence, stronger than the row-range inference above. The
[primary Jansen paper](https://research-portal.uu.nl/ws/files/63041030/Genome.pdf)
reports UKB `n=386533`, 28.3% prevalence, and linear or logistic regression
depending on phenotype. The raw OR and binary phenotype support a logistic
effect interpretation, but **the exact regression command and per-variant
case fractions remain undocumented**. Consequently, source `N_i` is
source-valid for a technical pilot, while full confirmatory admission and
bivariate promotion remain withheld pending model and overlap review.

### Count-field variability check (first 100,000 literal source rows each)

This is a source-only descriptive sample, **not** a whole-file statistic or a
pilot. Case fraction is calculated solely as `case_count/(case_count+control_count)`
for positive count fields.

| Trait | Minimum | Maximum | Mean | Distinct fractions rounded to four decimals |
|---|---:|---:|---:|---:|
| ADHD | 0.10954 | 0.41879 | 0.18120 | 486 |
| Parkinson | 0.04092 | 0.79772 | 0.35767 | 1,163 |
| Bipolar | 0.05539 | 0.44368 | 0.10608 | 1,408 |
| Schizophrenia | 0.36826 | 0.46523 | 0.40925 | 624 |

The Parkinson mean reflects the first 100,000 raw rows, including rows that
would later fail the frozen input QC; it is not a suggested LAVA fraction.
The PGC labels for bipolar and schizophrenia refer to effective counts, so
their ratios are descriptive of released fields, not actual variant case
fractions. This check rules out treating all four as having transparently
constant case/control composition. It does **not** determine any LAVA outcome.

## Prospective next actions (no result-based source selection)

1. **Insomnia first:** the official README binds `N_i` and OR to the exact
   release. Predeclare one source-faithful trait-only pilot arm using unchanged source
   OR/SE/P/effect alleles, literal `N_i`, and the UKB cohort case fraction
   `109402/386533 = 0.283034048839` as an explicit fixed-fraction
   approximation. The source lacks per-SNP case counts; document why that
   approximation is acceptable or hold admission. Freeze raw/harmonized
   hashes, case fraction, QC, allele mapping, representative loci independent
   of the 25 protected candidates, four partitions, gate, and acceptance
   criteria before materializing or reading outcomes. Keep LDSC N conventions
   separate from this LAVA-specific input. The pilot would assess technical
   and numeric feasibility only. Changed input requires matched pairwise LDSC
   overlap covariance before any bivariate LAVA.
2. **Parkinson:** request or find source-specific handling of proxy cases and
   whether published `N_cases`/`N_controls` are literal contributing analyzed
   participants for each meta-analytic statistic. Validate variable fractions,
   position/allele mapping, and overlap with UKB long sleep. The existing
   GP2 2025 public archive does not resolve the N/build issue without further
   metadata. Do not substitute count sums into LAVA yet.
3. **MDD:** pursue a source/method-supported balanced-equivalent PGC2025
   interpretation and new cohort overlap covariance. The previous 108/111
   pilot cannot choose the N model or justify a full trait run by itself.
4. **ADHD/bipolar/schizophrenia:** retain the exact sources, document cohort-wise
   meta-analysis and count semantics, and stage only if the four-trait rescue
   makes the family gate plausibly reachable. All source decisions precede
   association-level pilots.

### Method and provenance pointers

- `scripts/01_harmonize.py`, binary N branch: source `N_EFF`, then
  `N_EFF_HALF`, then `NCASE/NCONTROL`, then constant configured N_eff;
  literal `N` is used only for continuous traits. This explains the insomnia
  loss without altering historical LDSC intent.
- `config/gwas_schemas.tsv`, six original source rows, supplies exact header
  schema, build, allele direction, source-specific N conventions.
- `brain6/results/lava_binary_n_crosswalk_v1/binary_n_interpretation.md` and
  `brain6/results/lava_confirmatory_protocol_v1/sample_size_and_overlap_protocol.md`
  explain why effective N plus observed case fraction is a joint-model issue.
- `brain6/results/lava_multitrait_feasibility_v1/source_triage_20260926.md`
  and both terminal MDD pilot adjudications document MDD and GP2 limitations.
- The [Nalls primary study](https://pubmed.ncbi.nlm.nih.gov/31701892/)
  confirms proxy-case inclusion; the exact raw file has 10 count/statistic
  columns and no variant ID or INFO.

## Frozen native-N technical pilot outcome (v2)

The prospectively frozen and receipt-bound 88-locus, four-worker v2 technical
pilot is terminally held because one worker had a per-locus execution failure.
See [terminal pilot report](INSOMNIA_NATIVE_N_PILOT_V2_TERMINAL.md). Its
predeclared zero-failure advancement condition was not met, so it provides no
authorization for a full trait screen or confirmatory promotion. The direct
source N/OR/README evidence above remains valid, but no outcome changed the
source-admission class or the protected family gate.

The separately frozen [worker-2 postmortem](INSOMNIA_NATIVE_N_WORKER2_POSTMORTEM_V3_RESULT.md)
localized v2's failed cell to locus 950, which had zero rows in both the
canonical and native-N shard. This structural empty-input error has no bearing
on the direct source N semantics and does not change v2's failed decision.

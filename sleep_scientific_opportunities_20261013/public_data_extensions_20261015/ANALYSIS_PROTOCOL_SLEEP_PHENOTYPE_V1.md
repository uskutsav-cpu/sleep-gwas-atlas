# Frozen analysis protocol: sleep-duration measurement comparison

**Protocol ID:** `SLEEP-PHENOTYPE-RG-2026-10-15-v1`
**Status:** Frozen before acquisition or inspection of the GCST90824104 outcome file
**Primary source:** Portas et al. 2026 device-derived night-time sleep duration (GCST90824104)
**Scope:** One new three-source, source-aware sleep-phenotype comparison. This does not modify Atlas estimates, QC labels, thresholds, or historical testing families.

## Question and estimand

Among European-ancestry UK Biobank GWAS, is the genome-wide genetic correlation of deep-learning-derived, accelerometer-based night-time sleep duration with self-reported sleep duration different from its genetic correlation with the 2019 wrist-actigraphy sleep-period duration phenotype?

The single primary estimand is

`Delta = rg(GCST90824104, Dashti-2019 self-report) - rg(GCST90824104, Jones-2019 actigraphy)`.

The two-sided null is `Delta = 0`; the nominal type-I error is 0.05 for this one primary contrast. The estimand concerns shared common-variant genetic architecture between measurement-defined sleep phenotypes. It is not a disease association, causal effect, biological-mechanism claim, or independent replication test. The two secondary pairwise `rg` estimates and the self-report/legacy-actigraphy `rg` will be reported descriptively with uncertainty, without additional confirmatory hypotheses.

## Sources and admissibility

| Role | Exact source | Population and phenotype | Build and allele information | Independence and use |
|---|---|---|---|---|
| New device phenotype | GWAS Catalog GCST90824104, Portas et al. (2026), `GCST90824104.tsv.gz`; metadata MD5 `445dc180e3545392d171ba368a00e31a` | UK Biobank, European ancestry, combined sex, N=80,013; accelerometer-derived night-time sleep duration | GRCh37; unharmonized GWAS-SSF; allele orientation must be verified from literal data fields and harmonized against the frozen HM3 allele set | Public GWAS Catalog release; CC0 stated in release metadata. Participants overlap the other UK Biobank sources; not independent replication. |
| Subjective phenotype | Dashti et al. 2019, `sleepdur.sumstats.gz`; SHA-256 `7bba1d874e5c56842812ec3a973d5a14cc4779f4fefdb03a18bb4d19e6c2f4d2` | UK Biobank, European ancestry; continuous self-reported habitual sleep duration; N=446,118 | GRCh37; existing munged file has `SNP A1 A2 Z N` | Existing immutable project input. Cohort overlaps both device studies. |
| Legacy device phenotype | Jones et al. 2019, `accel_sleep_duration.sumstats.gz`; SHA-256 `724ba24f6294c9a3158735f084f7e4aa2c27d7afce38689e67c98e6c78aac958` | UK Biobank, European ancestry; rank-normalized wrist-actigraphy sleep-period duration; N=85,449 | GRCh37; existing munged file has `SNP A1 A2 Z N`; incomplete rows are excluded under the frozen complete-field rule | Existing immutable project input. Cohort overlaps the other two sources. |

The two existing munged files each contain 1,217,311 data rows; 1,179,618 self-report rows and 1,176,749 legacy-actigraphy rows have complete SNP, allele, Z, and N fields. Retain only rows passing the existing source QC and available in the frozen European HapMap3 (HM3) LDSC reference for all three traits. Require unique SNP IDs and valid, finite, positive N. Do not infer values for missing rows. Before fitting, record source-specific and three-way shared-SNP counts, allele conflicts, ambiguous alleles, and the number and proportion retained. If either existing input hash changes or its documented source identity cannot be reconciled, stop.

The three releases draw from UK Biobank. Exact person-level intersections are unavailable from these aggregate files. This is an explicitly overlapping-sample, within-UKB construct comparison. Estimate a cross-trait LDSC intercept for each bivariate fit using the standard LDSC behavior; report the intercept and its uncertainty. Do not describe any result as cohort-independent, external, or replicated.

## Frozen processing and estimator

1. Preserve raw downloads, the existing two munged inputs, and all original Atlas files. Write new source data and outputs only below the versioned `sleep_scientific_opportunities_20261013/public_data_extensions_20261015/` namespace or the existing SSD analysis namespace. Never overwrite an existing input or output.
2. Acquire only the one Catalog file, by the official release URL, to the SSD. Use resumable transfer from any verified partial checkpoint; validate final byte count, Catalog MD5, gzip CRC, and a locally computed SHA-256. Do not fetch the other three Portas traits.
3. Record the literal header and all relevant metadata before transforming the file. Map its GRCh37 variants to the existing HM3 reference by exact documented rsID or validated coordinate-and-allele key. Reject ambiguous strand alleles, duplicate/conflicting records, non-autosomal variants, invalid/non-finite effect or standard-error fields, and variants absent from HM3. Compute Z from beta/SE only when both are valid; retain the release's documented per-variant N where provided, otherwise use the documented study N=80,013 and label that assumption. Do not substitute a newer release or harmonized file.
4. Use the pinned LDSC implementation, European LD-score files, and weights already used by the project. Freeze their absolute paths and SHA-256 in the run receipt before estimation. Use the identical three-way SNP intersection, identical SNP order, one reference build, the standard LDSC regression settings, and exactly 200 delete-blocks for all three pairwise fits. Use standard allele alignment to the same reference effect allele. Do not use the MHC or variants excluded by the frozen project source-QC rules.
5. Run only the three pairwise bivariate LDSC estimates needed for this contrast. Retain full logs, point estimates, cross-trait intercepts, heritability and genetic-covariance deletion files, and all software/input/reference hashes. Do not launch unrelated Atlas or extension jobs.

## Paired uncertainty for the primary contrast

The three pairwise fits must use the same three-way SNP intersection and the same 200 genomic delete-block definitions. For each pair and each deletion block, reconstruct the deleted-block `rg` as `gencov_delete / sqrt(h2_trait1_delete * h2_trait2_delete)` from the LDSC output's aligned delete-block values. Confirm positive deleted-block heritabilities and exact block alignment; otherwise the primary contrast is not estimable under this protocol.

For the primary contrast, align the 200 `rg` delete-block values for device/self-report and device/legacy-actigraphy by block, subtract them blockwise, and calculate the delete-block jackknife standard error as `sqrt((B-1)/B * sum((Delta_b - mean(Delta_b))^2))`, with `B=200`. Report the point contrast from the full-data estimates, this paired standard error, the two-sided normal-approximation 95% interval and Wald P value. The paired deletion covariance is essential because both contrasts share the device GWAS and their variants; do not combine marginal standard errors as if the estimates were independent. Do not call covariance diagnostics calibrated until an independent numerical implementation reproduces the block alignment, contrast, SE, interval, and P value from the frozen raw deletion files.

The LD-block jackknife quantifies genomic sampling uncertainty conditional on the summary-statistic inputs and LDSC model. It does not identify participant-level covariance, eliminate population/sample-overlap bias, establish measurement reliability, or support causal direction. If the estimate is boundary-invalid, a deleted-block heritability is nonpositive, an intercept/QC gate fails under the existing project rules, variant/block alignment is incomplete, or independent arithmetic does not reproduce, report pairwise estimates descriptively and withhold inferential claims about `Delta`.

## Validation and interpretation

Independent numerical validation will use a separate, minimal calculation over the saved full-fit and delete-block files; it will not re-run the same LDSC pipeline as its only check. Validate hashes, parsing, block counts/order, the `rg` transformation, paired subtraction, jackknife formula, and reported interval/P value. Compare each full-fit pairwise `rg` to the output log. Preserve negative, null, unstable, or QC-failing outcomes.

Interpret `Delta` only as a difference in genome-wide common-variant sharing of operationally distinct sleep-duration phenotypes in overlapping UKB samples. Discuss differences in ascertainment, measurement model, sample size, phenotype transformations, and accelerometer-stage classification as plausible explanations, not demonstrated mechanisms. The Portas paper's independent PSG validation and previously published correlations are prior literature, not validation newly produced here. No locus-level, disease-level, ancestry-transferability, or causal follow-up is authorized by this protocol.

## Second comparison decision and stopping rule

The strongest independent-cohort candidate screened for this workstream is UKB insomnia symptoms versus FinnGen R13 registry insomnia (`F5_INSOMNIA`). Do not run it under the current source contract: the public R13 source is Finnish and GRCh38, its published summary-statistic schema lacks per-variant INFO and N, the prior h2 diagnostic used a qualified assumed effective N, and no authenticated Finnish genome-wide LDSC LD-score reference has been admitted. The official FinnGen LD browser is a local query service, not such a reference panel. The earlier frozen FinnGen protocol admits no pairwise `rg`. Reconsider only after an outcome-blind source gate establishes a licensed Finnish LD reference, variant-level information/N semantics, build/allele harmonization, adequate effective power, and a fresh protocol. A UKB-derived Portas source cannot repair independent sleep-cohort replication.

Proceed to estimation only after (a) this protocol is committed, (b) the official file is recovered or downloaded and passes its pinned integrity checks, (c) the exact old-input and HM3/LDSC reference hashes are recorded, (d) the existing process is not competing for memory/CPU, and (e) the native project's >=3 GiB internal-free-space and >=5 GiB SSD-free-space launch guard passes. If a gate fails, preserve the failure receipt and stop without substituting sources or thresholds.

## Reporting boundary

The only possible new scientific evidence from this experiment is the predeclared, calibrated within-UKB genetic-similarity contrast. A positive or negative result does not itself establish that the Atlas is more useful than existing resources, reveal a novel disease association, demonstrate biological mechanism, or support independent replication. No new discovery claim is allowed without additional, separately frozen and independently validated evidence.

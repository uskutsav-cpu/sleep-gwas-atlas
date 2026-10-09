# Frozen research decisions, v1

Frozen locally on 2026-10-08 before any newly acquired sleep-source disease correlations. This is technical research documentation, not manuscript text.

## New primary scientific question

**NO-GO at present. No novel primary hypothesis is admitted.** The independent novelty review ranks joint smoking/adiposity conditional genetic sharing as a future lead, but compatible independent data, calibrated covariance uncertainty, adequate power, and a distinct contribution have not passed their gates. No confirmatory biological experiment or local/molecular follow-up is authorized by this protocol. Significant historical associations are not a selection criterion for a new primary question.

## Native historical reproduction

Retain the exact atlas-v1.0 45-trait universe and 396-test family. Report the 372 historically primary rows separately from 24 sensitivity rows; 153 primary discoveries arise under the all-396 BH correction. The complete all-396 family contains 161 BH positives including eight sensitivity rows. Do not substitute the 372-only BH correction. Preserve the independent 100-phenotype/1,200-test extension family, its source versions, original h2/QC definitions, and BH correction. Reproduce the full 217 validation-candidate classification using the original 0.05/217 threshold, including sources that fail QC and candidates without sources. The 23 historical positives are qualified outcome-side validation because the original sleep inputs were reused.

Frozen tables were serialized from printed LDSC statistics. Native reruns retain unrounded return values and block-delete values, but a match to a rounded historical number cannot establish full-precision identity to an unavailable original return value. Compare rg/SE/intercepts at their recorded four-decimal precision, and P/Z at the exact historical printed precision. Keep native-full-precision BH and frozen-serialized BH as separate numerical audits.

## New independent clinical-insomnia transport validation

This is validation of an existing locked family, not the selection of a novel primary hypothesis. Source metadata are admitted conditionally before pair outcomes. The estimand is the global genetic correlation between **MVP clinically recorded insomnia (PheCode 327.4)** and each of the ten already registered FinnGen R13 outcomes among the 40 historical insomnia candidates. Clinical coding requiring at least two ICD instances is related to UKB field-1200 frequent self-reported complaints; these are not identical phenotypes. A positive result supports cross-phenotype transport, not exact replication of questionnaire insomnia or a causal mechanism.

The sleep source is the EUR-only GWAS Catalog accession **GCST90475826**, GRCh38, 78,566 cases/329,572 controls, N=408,138, effective N=253,768.6150468714. Exact compressed body: `GCST90475826.tsv.gz`, 575,504,178 bytes, upstream MD5 `2e9c12624b653aac527444fc426e037b`, ETag `224d7f32-630b35ab64241`. The full-body SHA-256 must be recorded before preprocessing or any pair outcomes. Source identity, primary methods, schema, access terms, and phenotype differences are in `reviews/provenance_v1_mvp_insomnia_priority.md` and the result-free 217-row admission table. Keep the raw file on this task's new SSD namespace and exclude it from Git.

Discovery sources contain UKB; the proposed validation sources are MVP EUR and FinnGen R13. Cohort compositions are distinct by study design. Exact individual cross-enrollment has not been measured and must remain unknown. The cohort-source distinction must not be turned into a claim of proven zero participant intersection. Reusing FinnGen R13 for the original FinnGen R9 sleep-apnea GWAS would not satisfy independence.

Preprocessing gates, prior to genetic-correlation results:

1. Match complete compressed-body bytes and MD5, pin SHA-256, and retain HTTP/source evidence.
2. Verify the literal effect allele and use beta=log(positive odds_ratio). Validate supplied SE against confidence-interval-derived log-OR SE and two-sided P; never assume its scale from the header. Use source-specific, documented per-variant effective N; examine missing/variable case-control counts before any constant-N fallback.
3. Apply the exact registered ancestry/build mapping, verified hg38-to-hg19 chain and HapMap3 identity/allele reference; exclude strand-ambiguous SNPs, invalid alleles, duplicate/ambiguous mappings, and invalid signed statistics. Inherit the historical replication preprocessing/QC rules rather than creating a more favorable threshold. Preserve absent per-variant INFO as a source limitation.
4. Require the historical h2 Z>=4 and intercept<=1.2 gates for sleep and outcomes. Recover or rerun the exact outcome inputs and verify their receipts. No newly estimated pair can be used to choose another outcome or release.

Power and uncertainty: h2-Z is an estimator-eligibility gate, not proof of power for rg. The prospective planning target is |rg|=0.15 with two-sided family alpha=0.05/217 and 80% power. Report the SE needed for this target and achieved detectable effects; a threshold-negative estimate without adequate power is not evidence of absence. No row is removed from the family on the basis of its observed rg, direction, or P.

Positive control: the native historical insomnia-BMI pair, under its original frozen source definition, is a computational reproduction control only. It does not validate clinical-phenotype equivalence. All ten eligible metadata candidates, including historically QC-failed outcomes, retain the same gates; no outcome is selected as an artificial biological negative control. No supported independent clinical negative-control phenotype has yet been admitted.

Success requires source/QC eligibility, the prespecified historical direction, and P<0.05/217=0.0002304147465437788. Classify successes only as **cohort-distinct clinical-insomnia transport validation** with the individual-overlap qualification. Every one of 217 rows stays in the ledger; ineligible/absent sources remain unavailable, not null genetic tests. If the sleep source fails h2/intercept or identity/schema gates, stop all ten pair tests. Preserve any failed or null eligible result. Do not change phenotypes, thresholds, source versions or methods after outcomes; any further proposal requires a new result-free version.

## Covariance work

Do not infer covariance from arbitrary correlation grids. Default 200-delete vectors from different pairwise SNP intersections are not automatically aligned genomic blocks. The original core 45-trait S/V matrices do not include the extension/FinnGen validation outcomes. GenomicSEM's fixed-diagonal V_Stand rescaling does not propagate uncertainty in estimated h2 denominators; it cannot be advertised as a complete ratio-rg sampling covariance. Any covariance correction must explicitly align genomic deletion boundaries, validate marginal SE against the native ratio jackknife, obtain independent implementation agreement, and keep the historical zero-covariance heterogeneity results as qualified diagnostics until that succeeds.

## Operational stopping rules

One native worker at a time; BLAS/OpenMP threads=1. Keep all large inputs on the SSD. Before launching each full native job, require at least 3 GiB free internally for swap headroom and 5 GiB free on the output volume. A failed resource preflight stops computation without changing scientific thresholds. Source download is limited to the single 575,504,178-byte admitted MVP file; no massive re-download or unbounded acquisition. Archive recovery uses one sequential reader and extracts only receipt-pinned inputs/QC/logs into a new directory.

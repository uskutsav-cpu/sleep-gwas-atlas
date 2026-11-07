# Whole-file checks of four physical-component GWAS

On 29 October 2026, four [Zenodo archives](https://zenodo.org/records/14011550) were downloaded to `/Volumes/Extreme SSD/brain6-public-component-gwas-20261029/`. All four completed decompression and matched their publisher MD5. The table compares source rows with uniquely identified rsIDs in the project's GRCh37 HapMap reference. Allele matches are unordered; mismatch counts are rows, not independent loci.

| Component | Total rows | Reference-ID rows | Coordinate matches | Unordered allele mismatches | Repeated reference IDs |
|---|---:|---:|---:|---:|---:|
| Low grip strength | 8,380,920 | 1,020,541 | 1,020,541 | 873 | 861 |
| Low physical activity | 8,380,973 | 1,020,540 | 1,020,540 | 873 | 861 |
| Slow walking speed | 8,380,994 | 1,020,540 | 1,020,540 | 874 | 862 |
| Weight loss | 8,381,014 | 1,020,537 | 1,020,537 | 874 | 862 |

The simple numeric checks detected no invalid rows. The roughly 7.36 million rows per file without a reference rsID were not assessed for coordinates or allele agreement in this comparison. Counts and hashes, including the exact reference hash, are in each component's `receipt.json`; scripts and per-chromosome tables are alongside them. These checks strengthen evidence of GRCh37 compatibility for the matched rows and file integrity, but do not establish effect direction, ancestry, covariates, exact phenotype coding, study overlap, or per-variant sample size. Repeated rsIDs and allele mismatches require allele-specific handling. No genetic associations were calculated and the inference hold remains in place. These separate GWAS cannot be treated as an exhaustion-deleted composite.

The [protocol](PROTOCOL.md) and earlier [full exhaustion audit](../zenodo_full_exhaustion_audit/RESULTS.md) describe the method and limitations. The four files are retained on the SSD; they are not bundled with this report.

## Follow-up on repeated IDs and mismatched alleles

A second full scan of each archive verified its MD5 and the SHA-256 from the first scan. For grip and activity, **862 of 873 mismatch rows** had a same-ID row with the reference allele pair; for walking and weight loss, **863 of 874** had one. Each file had **11 mismatch rows without a matching companion**. Every mismatched pair shared one allele with the reference pair. This pattern is compatible with multiple alleles at a site, but does not prove the representation or determine biological relevance. All repeated reference IDs were associated with a mismatch and an exact-pair companion.

The retained `duplicate_audit/allele_mismatches.tsv` and `receipt.json` in each component directory record the comparisons. Exact-pair matching avoids the specific rsID-only ambiguity identified here, but is only a candidate QC step: build, strand, palindromic variants and effect direction still need source-appropriate review. Do not guess the 11 unmatched pairs or promote these component files into cross-trait inference without the unresolved study metadata.

## Reported-statistic consistency

In a deterministic every-100th-row sample of each file (83,810 or 83,811 rows per component), the reported `LOG10P` agreed with the two-sided chi-square-one-degree-of-freedom tail computed from `CHISQ` within 0.000030 log10 units for every sampled row. No sampled row differed by more than 10% between `CHISQ` and `(BETA/SE)^2` using the denominator `max(CHISQ, 10⁻⁸)`; 4,018–6,112 rows per file differed by more than 0.1%. Small differences may reflect rounding or a different test; these comparisons do not identify the method or establish calibration. The sample is ordered and variants are correlated, so it is not a random sample or independent test set. Source SHA-256 was rechecked, and the [script](statistic_check.py) and component-specific `statistic_check.json` files record counts. Effect-allele coding, phenotype and cohort metadata remain unverified.

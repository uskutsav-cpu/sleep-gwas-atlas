# Fried Frailty Score full GWAS acquisition and QC audit

**Acquired:** 2026-09-26 UTC
**Status:** Full summary-statistics file acquired and checksum-verified; coordinate build is now empirically supported as GRCh37/hg19. The source-consistent A1/BETA convention is documented as a FUMA/BOLT inference; exact participant overlap and other analysis gates remain open.

## Source and transfer

Ye et al. (2023; PMID 36928559; DOI 10.1007/s11357-023-00771-z) describe the UK Biobank Fried Frailty Score (FFS) GWAS in 386,565 European-ancestry participants. The score is ordinal from 0–5, and the paper reports BOLT-LMM v2.3.2 single-variant testing with age, sex, and 20 genetic principal components as covariates. The article-linked Figshare share page identifies `Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv`, displays CC BY 4.0, and reports a 609.63 MB file.

The file was retrieved using Figshare's documented `ndownloader.figshare.com/files/{file_id}?private_link=...` route. A byte-range response reported 639,246,335 total bytes. The complete transfer was written to a separate partial filename, verified, then atomically promoted to the registered target.

## Integrity and stream validation

- Registered path: `frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv`
- Bytes: 639,246,335 (matches the published range total)
- SHA-256: `912e290b2064999c7147cbe3597f87cbcb14b25d387a8b6db2f25862e7274e55`
- MD5: `1c6a2f40537f68c2d06d0c55a4ce1dee` (matches the S3 ETag returned through the Figshare downloader)
- Header: `CHROM POS A1 A2 BETA P SNP SE OR`
- Data rows: 8,883,488
- Chromosomes: 1–22; all positions positive
- Stream scan found no wrong-width rows, invalid numeric values, non-finite effects/statistics, out-of-range P values, non-positive SE/OR, or noncanonical alleles.
- External bounded-memory sort found zero duplicate `(CHROM, POS, A1, A2)` keys.
- For all rows, the file's `OR` column equals `exp(BETA)` within the explicit numerical tolerance in the scan. It is therefore not treated as an independently estimated odds ratio.

## Eligibility decision and remaining source checks

The downloaded file is registered in `manifests/all_acquired_resources.tsv` and `manifests/file_checksums.tsv`. The streaming coordinate audit in `analysis/fried_ffs_grch37_coordinate_concordance_2026-09-26.md` found 1,162,260 shared coordinates with the pinned GRCh37/HapMap3 map, including 1,162,055 matching rsIDs and zero sort violations. This is strong exact-file empirical support for GRCh37/hg19. The source-consistent signed-effect interpretation is now documented in `analysis/fried_ffs_effect_coding_recheck_2026-09-26.md`: infer `A1` as the effect allele from the exact FUMA-labelled file/header and official FUMA/BOLT conventions. This remains an inference rather than author-sidecar confirmation. The audit does not establish per-variant sample sizes or exact participant intersections. Do not claim independent replication; overlap-dependent comparisons remain qualified or blocked under the frozen plan.

The study is UK Biobank-based, so cohort overlap is expected for the 11 UKB-based sleep sources; exact intersections remain unknown. It must remain a distinct construct from the primary Frailty Index and cannot be described as independent replication. No LDSC, LAVA, or other inferential analysis was run on this file. The frozen plan and thresholds were unchanged.

## References

- Ye et al. 2023, [primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10651618/).
- [Article-linked Figshare share page](https://figshare.com/s/6683396c68807fe4e729), file 38980967.
- Figshare [official API documentation for private-linked downloads](https://docs.figshare.com/old_docs/api/articles/).

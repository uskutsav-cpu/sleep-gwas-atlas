# Fried FFS coordinate-build concordance audit — 2026-09-26

## Question and decision

Does the exact article-linked Ye et al. full-statistics file use GRCh37/hg19 coordinates? Yes: a streaming coordinate comparison against the project's pinned GRCh37/HapMap3 variant map provides strong exact-file support for GRCh37/hg19. This resolves the registered coordinate-build field. It does not resolve effect-allele direction, per-variant sample size, or UK Biobank participant overlap; the file remains ineligible for harmonization and inference until those gates are addressed.

## Reproducible comparison

The auditor is `frailty_paper/scripts/56_audit_fried_build_concordance.py`. It filters the summary file to rsID-labelled single-nucleotide variants, groups any rows sharing a chromosome and 1-based position, then merge-joins the sorted coordinate groups. Within each shared coordinate group it pairs exact rsIDs and compares unordered allele pairs. It stores only the current coordinate groups in memory. Both streams had zero sort violations.

- Summary file: `frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv`; SHA-256 `912e290b2064999c7147cbe3597f87cbcb14b25d387a8b6db2f25862e7274e55`.
- Reference map: `/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/hm3_grch37_variant_map.tsv.gz`; SHA-256 `6775a7a0d3ca90dc74e472180b1d77103bc238129c4f969358f46307e5c306b4`. Its adjacent provenance identifies GRCh37/hg19 and records its source-file hashes; provenance SHA-256 `5f786febaba44c8f561c2c1b4402901c91fa7666e03fc3bf68fe73584d4c0516`.
- Compared summary variants: 7,981,276 rsID-labelled biallelic A/C/G/T rows.
- Shared coordinate groups: 1,162,260; exact coordinate-and-rsID variant matches: 1,162,055; unmatched variant rows within shared coordinate groups: 872 (including positions with multiple rows).
- Of same-coordinate/same-rsID variants, 351,435 used the same allele order and 810,404 the reversed order; 216 had allele-pair mismatches. Allele labels in the reference map were used only for unordered pair comparison, not to infer effect direction.
- Machine-readable results: `analysis/fried_ffs_grch37_coordinate_concordance_2026-09-26.json`.

The high concordance across over 1.16 million common variants, together with study-level GRCh37 metadata already identified for GWAS accession GCST90295968, supports GRCh37/hg19 for this exact file. This is an empirical concordance result against the locked map; it is not a liftover or proof that every row is valid.

## Effect coding remains unresolved

Ye et al. report BOLT-LMM association testing for this GWAS, and the official BOLT-LMM manual specifies that its standard `ALLELE1` field is the effect allele for `BETA`. The acquired file instead has columns `A1`, `A2`, and `BETA`, without a source-specific statement mapping `A1` to BOLT's `ALLELE1`. Therefore the method documentation is supportive but not sufficient to assert the exact file's signed-effect convention. No effect signs were changed, and no harmonization or downstream inference was performed. The `OR` field remains a derived `exp(BETA)` value per the acquisition audit, not an independent odds-ratio estimate.

## Sources

- Ye et al. primary article, which reports the UK Biobank FFS GWAS, BOLT-LMM analysis, and links the full summary statistics: https://pmc.ncbi.nlm.nih.gov/articles/PMC10651618/
- Official BOLT-LMM manual, which defines the standard `ALLELE1` effect-allele convention: https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html
- Official PGS Catalog source record for the corresponding study accession, with study-level GRCh37 metadata: https://www.pgscatalog.org/score/PGS005229/

# HFRS publisher-supplement acquisition audit

Audit timestamp: 2026-09-24 03:12 UTC.

## Source and acquisition

The primary article is Mak et al., *Nature Aging* 5:1589–1600 (2025), DOI [10.1038/s43587-025-00925-y](https://doi.org/10.1038/s43587-025-00925-y). Its public article page links the workbook [43587_2025_925_MOESM1_ESM.xlsx](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs43587-025-00925-y/MediaObjects/43587_2025_925_MOESM1_ESM.xlsx) as Supplementary Tables 1–19 and states the article is CC BY 4.0. Retrieved over ordinary HTTPS from the publisher on 2026-09-24 and archived unchanged at `frailty_paper/data/gwas/hfrs_mak_2025/mak_2025_hfrs_supplementary_tables.xlsx`. The 1,059,282-byte XLSX passed ZIP integrity (`ZipFile.testzip()` returned no corrupt member); SHA-256 is `123c340cf063e2ef618b8ef543dedea6b7c8ef69008488f01027eedbffc04625`.

## Table scope verified from the workbook

The workbook contains 20 sheets: Contents and ST1–ST19. The Contents sheet titles Supplementary Table 1 as variants associated with HFRS at P<5×10⁻⁸ in FinnGen with UK Biobank replication and METAL meta-analysis; ST2 gives the corresponding HFRS-without-dementia results. Read-only inspection with openpyxl 3.1.5 found 29 columns in both result sheets, including FinnGen beta/SE/P, UKB beta/SE/P, and METAL Z/P and direction. ST1 contains 1,588 nonempty result rows; ST2 contains 492. These are thresholded reported-hit rows, not all variants or genome-wide complete summary statistics. Their presence does not make LDSC/LAVA possible and does not establish that the full custom-score summary files are public.

The article reports FinnGen R12 discovery N=500,737 and UK Biobank replication N=407,463. The publisher workbook does not state a genome build in its table headers or introductory metadata; do not infer the build from coordinates. UKB is a known cohort-level overlap for downstream sleep/FI comparisons; exact participant intersection is unresolved. The article’s Data Availability section refers generally to post-embargo FinnGen summary-statistic releases but names no HFRS-specific accession or direct complete file. FinnGen’s [current results page](https://www.finngen.fi/en/access_results) describes DF13 and says bulk summaries follow an online form and emailed instructions; the standard DF13 endpoints are not evidence that the custom R12 HFRS measures are available. No form was submitted and no access was bypassed.

## Use boundary

Register the source as publisher supplementary GWAS hit tables, not as full GWAS summary statistics. Use only for qualified reported-locus cross-reference after build, variant, allele and source checks. Do not use it for genome-wide correlation, local-correlation scans, or independent UKB replication claims. The exact full HFRS and no-dementia files remain an open source/access gate.

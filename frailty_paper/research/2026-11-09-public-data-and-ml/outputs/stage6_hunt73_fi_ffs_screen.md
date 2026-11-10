# Exploratory circadian variant screen in FI and Fried GWAS

We scanned two existing public GWAS-derived files for the 73 rsIDs from the HUNT3 circadian-variant insomnia panel. The FI file is the project's harmonized, filtered Atkins UK Biobank/TwinGene summary file; the FFS file is Ye et al.'s UK Biobank Fried-score GWAS. This is a post hoc panel screen of marginal frailty associations, **not validation of a score, a sleep–frailty joint test, or a colocalization analysis**. Both frailty discovery sources involve UK Biobank and the FI includes sleep-related items.

| Source | HUNT rsIDs found | Same unordered allele letters | Tests with numeric p |
|---|---:|---:|---:|
| FI, filtered harmonized file | 68/73 | 67 | 67 |
| FFS | 72/73 | 71 | 71 |

For this descriptive set of 138 p-values pooled across both outcomes, three pass Benjamini–Hochberg 5%: FI rs2269457 (HUNT gene label NR1D1; p=0.0002728; adjusted p≈0.0228), FI rs2585408 (PER1; p=0.0003311; adjusted p≈0.0228), and FFS rs10462028 (CLOCK; p=0.000520; adjusted p≈0.0239). These are frailty associations in a postselected clock-gene panel. The SNPs and phenotypes may be correlated, and this adjustment does not include the seven HUNT sleep phenotypes per variant or the wider project's searches. The labels are the HUNT supplement's labels, not proven causal genes.

For context, the HUNT workbook's seven *unadjusted* sleep-subtype p-values for these three variants range 0.0735–0.9517, 0.0190–0.9983, and 0.2205–0.9447 respectively. The nominal 0.019 is one of seven selected contrasts for one selected variant, with no study-wide correction calculated here. These numbers do **not** establish shared sleep–frailty risk. We did not compare effect directions: HUNT's counted-allele semantics, the original FI source, the FFS convention, exact genome build/position and multiallelic context require source-specific validation. An rsID and two matching allele letters alone cannot settle those questions. A missing or mismatched row is not a null association.

The [row-level screen](stage6_hunt73_fi_ffs_allele_screen.csv), [source checksums and counts](stage6_hunt73_fi_ffs_receipt.json), and [scan code](../work/hunt73_fi_ffs_screen.py) document the result. We corrected an initial receipt-counting error caused by defaultdict key creation before writing this review; the table's row-level coverage was unaffected. Reproduction requires the cached source files at the paths recorded in the script. No participant-level data were used. The [FinnGen HFRS panel screen](stage6_hunt73_hfrs_circadian_screen.md) evaluates a different hospital-based frailty phenotype. Neither screen validates the 15- or 22-variant score.

## GRCh37 reference check for the three highlighted variants

On 29 October 2026, the [Ensembl GRCh37 variation API](https://grch37.rest.ensembl.org/) returned one GRCh37 mapping for each highlighted rsID. The FI/FFS summary position and HUNT allele letters agree with the API mappings:

| Variant | GRCh37 position | Ensembl allele string | Source summary position and pair |
|---|---|---|---|
| rs2269457 | 17:38254689 | T/C | FI 17:38254689, T/C |
| rs2585408 | 17:8043753 | C/T | FI 17:8043753, T/C |
| rs10462028 | 4:56298300 | G/A | FFS 4:56298300, G/A |

The HUNT supplement does not itself list these coordinates, and the API mapping does not verify HUNT's effect-allele semantics or the source cohorts' variant handling. These are single-base, two-allele strings in the returned Ensembl mapping, not proof of absence of other alleles in every cohort. This improves reference consistency but does not verify signed effects, establish colocalization, establish an insomnia association or validate a score. Raw API responses are `work/hunt_rs2269457_ensembl37.json`, `work/hunt_rs2585408_ensembl37.json`, and `work/hunt_rs10462028_ensembl37.json`.

An additional [FinnGen R13 insomnia lookup](stage6_hunt3_finngen_insomnia_review.md) matched all three highlighted variants by GRCh38 position, rsID and reference alleles. Insomnia association estimates for all three had confidence intervals crossing one; this is a postselected follow-up and does not exclude small effects or establish shared risk.

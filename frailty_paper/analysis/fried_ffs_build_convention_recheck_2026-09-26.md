# Fried Frailty Score build and allele-convention recheck — 2026-09-26

## Question

Can the exact Figshare summary-statistics object acquired for Ye et al. (2023), `Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv`, now be treated as GRCh37 with confirmed effect-allele semantics?

## Evidence checked

- The primary Ye et al. article describes the phenotype as ordinal FFS 0–5 in 386,565 European-ancestry UK Biobank participants. It reports single-variant GWAS using BOLT-LMM v2.3.2 and says the full statistics are at the linked Figshare share. It does not, in the accessible methods text, identify the genome build or define the acquired file's A1/A2 direction.
- A 2024 cross-trait paper identifies the Ye et al. (2023) frailty phenotype GWAS as FP, reports that the included GWAS summary statistics are GRCh37, and states it used LDSC `munge_sumstats.py` with signed Z statistics. Its text does not name the exact Figshare share or the acquired filename, nor explain the Ye file's allele coding.
- The acquired object has columns `CHROM POS A1 A2 BETA P SNP SE OR`. In the byte-validated file, `OR == exp(BETA)` for every row. That arithmetic identity is not source documentation and does not establish that OR is an independently estimated odds ratio or resolve the BETA sign's allele.

## Decision

Retain the exact-file genome build as UNKNOWN for harmonization. GRCh37 is corroborated at the study level by the later cross-trait paper, but the chain from that statement to this exact full-statistics Figshare object is not explicit enough for a locked source-specific build declaration. Retain effect-allele semantics as unresolved; neither the column names nor the BETA/OR relationship substitutes for the source's coding key. Do not harmonize, flip alleles, or analyze this file pending a source statement or an exact-object metadata record that defines both coordinate build and effect allele.

The acquisition and structural QC remain valid. The exact UK Biobank participant intersections with the sleep sources also remain unknown.

## Sources

- Ye et al. (2023), primary article, PMID 36928559 / DOI 10.1007/s11357-023-00771-z: [PubMed record](https://pubmed.ncbi.nlm.nih.gov/36928559/) and [PMC article](https://pmc.ncbi.nlm.nih.gov/articles/PMC10651618/).
- El-Saadi et al. (2024), cross-trait meta-analyses: [PMC article](https://pmc.ncbi.nlm.nih.gov/articles/PMC11069310/), section “GWAS summary statistics.”
- Article-linked exact-statistics share: [Figshare share](https://figshare.com/s/6683396c68807fe4e729) (web open returned HTTP 403 in this check; local downloaded bytes remain checksum-verified).

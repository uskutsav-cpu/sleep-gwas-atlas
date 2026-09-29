# Fried Frailty Score GWAS source and access audit

Audit date: 2026-09-23. Purpose: identify the authoritative full-summary-statistics source for the locked Fried Frailty Score (FFS) endpoint, verify its study/phenotype metadata, and determine whether acquisition is currently possible through the cited route.

## Exact study and phenotype

The primary study is Ye et al. (2023), *A genome-wide association study of frailty identifies significant genetic correlation with neuropsychiatric, cardiovascular, and inflammation pathways* (GeroScience; PMID 36928559; DOI [10.1007/s11357-023-00771-z](https://doi.org/10.1007/s11357-023-00771-z)). The authors analyzed the Fried Frailty Score as an ordinal score from 0–5, calculated from five criteria, in 386,565 unrelated participants of genetically confirmed European ancestry from UK Biobank. Their methods report exclusion of 35,588 participants for missing one or more criteria; BOLT-LMM v2.3.2; covariate adjustment for age, sex and 20 genetic PCs; and variant filters of imputation INFO >0.9, HWE P>1e-6, MAF>0.01 and genotyping missingness <0.1. The paper reports an independent HRS validation cohort of 9,720 participants, but the linked Figshare full-statistics item is specifically the UKB FFS file, not an independent sleep–frailty result.

The paper's Data Availability section links to a Figshare share item. Its visible item metadata identifies the file as `Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv`, describes it as GWAS summary statistics of Fried Frailty Score in UK Biobank, gives a size of 609.63 MB and displays a CC BY 4.0 license. A browser preview displayed columns `CHROM`, `POS`, `A1`, `A2`, `BETA`, `P`, `SNP`, `SE`, and `OR`. This is consistent with the intended FFS discovery resource, but the preview is not the full dataset and does not independently establish genome build or all source/QC metadata.

## Access status

The official article's cited Figshare share page loads and exposes a limited row preview and file metadata. The prior direct file request for Figshare file 38980967 returned HTTP 202 with `x-amzn-waf-action: challenge` and no file body. A standard browser download attempt in this audit did not yield a confirmed completed download within the browser operation window; no local copy was verified or promoted. The existing zero-byte local target remains an acquisition failure, not an acquired GWAS. No CAPTCHA was solved, no challenge bypassed, and no credentials or forms were used.

## Eligibility and remaining checks

The exact phenotype is a strong candidate for the locked Fried construct-replication family because it uses a five-component Fried score and matches the reported sample size in the Catalog metadata. Keep it unavailable until the full file can be obtained through a normal, verified route. On acquisition, record its byte count, SHA-256, decompression/readability, complete schema, chromosome/position build, allele/effect semantics, sample-size fields and source checksum before considering harmonization. UK Biobank overlap with most sleep GWAS is expected; exact participant overlap remains unknown. The paper's HRS validation is phenotype-side validation only and does not supply a paired independent sleep–FFS association.

## Sources

- Ye et al. full text and Data Availability statement: [PMC article](https://pmc.ncbi.nlm.nih.gov/articles/PMC10651618/)
- Figshare share item linked by the article: [Figshare item](https://figshare.com/s/6683396c68807fe4e729)
- GWAS Catalog metadata for GCST90295968 records 386,565 European discovery and 9,720 European replication participants, with `fullPvalueSet: false`.

## GWAS Catalog and PGS Catalog cross-check (2026-09-24)

The live GWAS Catalog page for `GCST90295968` explicitly reports “Full Summary Statistics: Not available.” The PGS Catalog separately exposes `PGS005229` (`FFS_Frailty`), a GRCh37 SBayesR score with 2,761,868 beta-weighted variants whose source GWAS is the 386,565-person European UK Biobank Fried Frailty Score study. Foote et al. (2025) report evaluation of this FFS-derived score against Frailty Index outcomes in independent older-adult cohorts: ELSA (`n=7,181`), PISA (`n=3,265`), and LBC1936 (`n=1,005`). The PGS Catalog reports, for example, an ELSA FFS-score association with FI (OR 1.08; beta 0.076; R² 0.0138; FDR-corrected p=2.48×10⁻⁶). Sources: [GWAS Catalog GCST90295968](https://www.ebi.ac.uk/gwas/studies/GCST90295968), [PGS Catalog PGS005229](https://www.pgscatalog.org/score/PGS005229/), and [Foote et al., Nature Genetics (2025)](https://doi.org/10.1038/s41588-025-02269-0).

This is cross-cohort polygenic prediction of frailty-related outcomes, not a full FFS GWAS release and not replication of any sleep–frailty association. The score file contains shrunken prediction weights rather than the full association estimates, standard errors, and P-values needed for the planned LDSC analysis. It must not be substituted for the inaccessible Ye et al. discovery summary statistics or treated as pairwise sleep–FFS replication. The Ye discovery and the PGS source cohort are UK Biobank; exact overlap of the PGS evaluation samples with any sleep resources is not inferred.


## Third ordinary browser download recheck (2026-09-23)

The visible Figshare `Download (609.63 MB)` link was activated through the in-app browser. A local file appeared in the standard Downloads directory as `Unconfirmed 591651.crdownload`; at inspection it was 6,296,188 bytes and its first bytes were consistent with the advertised tab-delimited GWAS columns and rows. The file remained the `.crdownload` partial and unchanged at the same size across more than 35 seconds of observation. The project-local FFS placeholder remained zero bytes. No complete file, expected byte count, final filename, checksum, or full-file validation was obtained.

This is evidence that the ordinary link began a transfer, not evidence of successful acquisition. The incomplete Downloads item is left untouched because no authoritative completed/canceled transfer state was available; it must not be used as an analysis input. No challenge, CAPTCHA, login, or access control was bypassed. The source remains unavailable for harmonization and analysis pending a complete normal-route transfer or author-provided copy.

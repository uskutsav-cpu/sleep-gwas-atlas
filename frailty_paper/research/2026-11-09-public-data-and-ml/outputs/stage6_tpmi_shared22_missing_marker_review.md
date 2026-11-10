# TPMI missing-marker follow-up for the 22-marker score

On 28 October 2026 the five 22-marker variants without an exact TPMI rsID autocomplete match also produced no matches when searched by their reference-checked GRCh38 positions. A positive-control coordinate search for rs12536335 did return the expected variant. Search behavior is not a comprehensive query of the complete summary file, so this does not establish absence from the target GWAS.

I checked the corresponding GRCh37 variant records against the public 1000 Genomes phase 3 East Asian (EAS) reference panel. Each had one exact biallelic position, rsID and allele-pair match, with 504 called individuals (1,008 chromosomes). The results count the allele specified in the 22-marker score:

| Variant | Counted allele | Counted copies | Minor allele copies |
|---|---|---:|---:|
| rs34290943 | T | 1,006 | 2 |
| rs9536347 | T | 0 | 0 |
| rs9845387 | C | 1,008 | 0 |
| rs12524473 | A | 1,006 | 2 |
| rs17312861 | G | 1,007 | 1 |

All five have a panel minor-allele frequency below 1%. [Chen et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12675292/) report a post-imputation frequency filter for their TPMI analysis. Low frequency is therefore a plausible explanation for missing browser results and a warning about transport of the particular score to this population. The EAS reference panel is not the TPMI cohort: these observations do not establish target frequencies, imputation quality, filtering, variant absence, or score performance. Two variants show no variation in this sample, which is not proof that they are monomorphic in TPMI or all East Asian populations. Multiallelic records in other resources and differences in ancestry composition warrant further checks.

The [machine-readable receipt](stage6_tpmi_shared22_missing_eas_frequency.json) and [script](stage6_tpmi_shared22_missing_eas_frequency.py) record the inputs and counts. The remote VCF was accessed by indexed region and was not independently full-file checksum verified. Source coding and the missing full-precision associations remain unresolved; neither a score-level association nor individual-level validation follows from these reference frequencies.

TPMI requests acknowledgement and citation of [Chen et al.](https://doi.org/10.1038/s41586-025-09350-y) and [Yang et al.](https://doi.org/10.1038/s41586-025-09680-x).

# Reference frequencies for three unverified TPMI markers

I checked the three unresolved 15-marker prototype sites in the 1000 Genomes phase 3 East Asian (EAS) reference panel. All three had one biallelic rsID/position/allele-pair match in the hg19 reference and 504 called individuals:

| Variant | Project counted allele | Counted-allele copies / 1,008 | Counted-allele frequency |
|---|---|---:|---:|
| rs113851554 | T | 0 | 0 |
| rs78206187 | G | 1 | 0.0010 |
| rs323509 | A | 69 | 0.0685 |

The TPMI paper reports excluding variants with MAF ≤ 0.01 and INFO ≤ 0.7 during post-imputation QC ([Chen et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12675292/)). The two very low EAS-panel frequencies offer a **possible** explanation for the first two search misses, but they are not TPMI allele frequencies, and the panel includes several East Asian populations rather than the target cohort alone. No conclusion about TPMI inclusion, imputation quality or absence can be made without its actual data. The rs323509 frequency is above that threshold in this panel, so low reference frequency does not explain its missing browser result; alternate representation and target QC remain possible and unverified. For rs323509 and rs113851554, other sources list additional alleles; the calculations apply only to the biallelic representations in this reference VCF.

These reference frequencies do not provide an outcome association or validate any score. Reproduction is in `stage6_tpmi_missing_eas_frequency.py` and its JSON output; the remote reference VCF was not independently full-file checksum verified.

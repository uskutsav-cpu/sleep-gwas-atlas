# TPMI browser effect-allele comparison with the published table

The authors' [Supplementary Table S4](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-025-09350-y/MediaObjects/41586_2025_9350_MOESM3_ESM.xlsx) labels non-effect allele (NEA), effect allele (EA), effect-allele frequency, beta, SE and p. For three selected records, I compared those fields with the publicly displayed TPMI PheWeb variant page on 28 October 2026. Browser values are rounded.

| Phenotype and marker | Table NEA / EA, beta (SE), p | Browser REF / ALT, effect (SE), p |
|---|---|---|
| Organic or persistent insomnia 327.41; [rs78324452](https://pheweb.ibms.sinica.edu.tw/variant/2-171723595-C-T) | C / T, +0.0979 (0.0178), 1.89e−8 | C / T, +0.098 (0.018), 1.9e−8 |
| Sleep-related movement disorders 327.7; [rs2814898](https://pheweb.ibms.sinica.edu.tw/variant/6-38420930-A-G) | A / G, −0.2102 (0.0323), 4.03e−11 | A / G, −0.21 (0.032), 4.0e−11 |
| Viral hepatitis 070; [rs74124815](https://pheweb.ibms.sinica.edu.tw/variant/1-160453302-G-A) | G / A, −0.0977 (0.0131), 4.59e−14 | G / A, −0.098 (0.013), 4.6e−14 |

All three published effect alleles are the browser's displayed alternate allele, and signed effects, SEs, p-values and case/control counts agree to browser rounding. This is direct source-specific evidence supporting an ALT-effect interpretation for these examples, including one insomnia phenotype and both positive and negative signs. The examples were selected from the paper's significant signals and are not project score variants. They do not trace the complete production conversion, guarantee every phenotype or variant was transformed consistently, resolve multiallelic records, or confirm counted-allele orientation of the project weights. The target score variants should continue to be marked unverified until their exact source and transformation are checked. The published insomnia signal does not pass the authors' across-phenotype Bonferroni threshold in S4 and is not independent evidence of a shared insomnia–frailty causal variant.

The [supplement extraction receipt](stage6_tpmi_supplement_insomnia_audit.json) records the workbook hash; S4's three source records can be extracted by exact phenotype and rsID. This comparison does not provide individual-level validation, target LD, or a score effect. TPMI requests acknowledgement and citation of [Chen et al.](https://doi.org/10.1038/s41586-025-09350-y) and [Yang et al.](https://doi.org/10.1038/s41586-025-09680-x).

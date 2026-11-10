# Public general-frailty score as a future benchmark

[PGS005221](https://www.pgscatalog.org/score/PGS005221/) is a public 943,126-row, GRCh37 genome-wide score from [Foote et al.](https://doi.org/10.1038/s41588-025-02269-0). It uses SBayesR weights for a latent general frailty factor derived using genomic structural equation modeling across 30 frailty-related traits. It is a candidate **frailty comparator** if a suitable independent participant dataset becomes available, not a substitute for the project's FI and Fried-specific scores.

In the source study, the general-factor score was evaluated against FI in ELSA (n=7,181), PISA (n=3,265) and LBC1936 (n=1,005). **Correction, 29 October 2026:** I previously presented the PGS Catalog's variance-explained values without checking the corresponding rows of the paper's supplement. The [paper's Supplementary Table 44](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-025-02269-0/MediaObjects/41588_2025_2269_MOESM4_ESM.xlsx) and [PGS Catalog](https://www.pgscatalog.org/score/PGS005221/) do not fully agree:

| Cohort, FI outcome | Catalog R² | Supplement Table 44 GF-PRS R² |
|---|---:|---:|
| ELSA | 0.009 | 0.01 |
| PISA | 0.006 | 0.02 |
| LBC1936 (Wave 1 in supplement) | 0.0002 | 0.01 |

The ELSA values are consistent with rounding; the PISA and LBC values are not. The source descriptions checked do not explain the differences, which could reflect different specifications or an error. Do not use those discrepant values for a numerical ranking until reconciled. The paper's [figure and caption](https://genepi.qimr.edu.au/contents/p/staff/Foote%20et%20al.,%20Uncovering%20the%20multivariate%20genetic%20architecture%20of%20frailty%20with%20genomic%20structural%20equation%20modeling.pdf) describe the cohort comparisons, which concern the published score, not our prototypes. No participant-level data were obtained here.

## Project candidate overlap

I downloaded the complete scoring file and compared genomic position and unordered allele pair with the project's GRCh37 candidate records. The script verified 943,126 rows. There were **0/15** and **3/22** direct position-plus-pair matches; the three latter matches also share the reported rsID. The [candidate rows](stage6_frailty_general_pgs_matches.csv), [script](stage6_frailty_general_pgs_comparison.py), and [receipt](stage6_frailty_general_pgs_receipt.json) document the check. The candidate coordinates are inferred or assumed in the input files, and matching a position and allele pair does not settle multiallelic contrasts. Nonmatching rows can still be in LD with score variants; low exact overlap does not mean little genetic overlap, and it does not evaluate the project scores.

The general-factor model **includes insomnia among its input traits** and uses UK Biobank data also represented in project discovery. An association of this score with insomnia would therefore not be independent evidence of sleep–frailty sharing. Its score-level associations also cannot be assigned to the few matched markers. The paper's cross-cohort analyses used genotype data; public weights and published outcomes do not give us participant-linked data. Before a head-to-head comparison, verify licensing, independent cohort/relative overlap, allele coding and coverage, target-outcome definitions, and test prespecified models in the same participants. A sleep-excluded frailty comparator, if available, would help address phenotype overlap.

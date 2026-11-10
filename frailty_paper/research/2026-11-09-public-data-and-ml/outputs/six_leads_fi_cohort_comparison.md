# FI cohort consistency for six cross-trait candidates

The FI meta-analysis hides potentially important cohort differences for several of the six candidate variants. In particular, **rs28528997**, the MAD1L1-region candidate with a strong FFS association, has opposite estimated effects in the original UK Biobank and TwinGene FI components under the project’s inferred TwinGene sign adapter. The same is true for rs12536335 and rs2815764. These are exploratory cohort comparisons; they do not establish biological differences, or a causal effect of sleep on frailty.

All effects below are aligned to the allele shown and are on the separately quantile-normalized FI scales. Differences are TwinGene minus UK Biobank, with approximate 95% Wald intervals under an assumption of independent cohort estimates. These p-values are **unadjusted** and were calculated after candidate selection and prior inspection of the data.

| Variant; counted allele | UKB beta (SE) | TwinGene beta (SE), inferred adapter | Difference (95% CI) | Difference p |
|---|---:|---:|---:|---:|
| rs28528997; G | −0.0213 (0.00436) | +0.0549 (0.0170) | +0.0763 (+0.0419 to +0.1107) | 1.36 × 10⁻⁵ |
| rs12536335; G | −0.0189 (0.00345) | +0.0392 (0.0145) | +0.0581 (+0.0288 to +0.0874) | 0.000101 |
| rs2815764; G | +0.0187 (0.00430) | −0.0342 (0.0183) | −0.0529 (−0.0898 to −0.0161) | 0.00491 |
| rs1909650; A | −0.0181 (0.00350) | +0.00285 (0.0143) | +0.0210 (−0.00786 to +0.0498) | 0.154 |
| rs817509; G | +0.0135 (0.00538) | −0.0131 (0.0212) | −0.0266 (−0.0695 to +0.0163) | 0.224 |
| rs2155292; G | +0.0191 (0.00348) | +0.0132 (0.0145) | −0.00592 (−0.0351 to +0.0232) | 0.691 |

For rs817509, which showed aligned nominal signals across six trait sources, the TwinGene point estimate has the opposite sign from UKB but is imprecise; the difference interval includes zero. That cannot establish either consistency or inconsistency of true effects. For rs2815764 the TwinGene estimate alone has p=0.062, yet the direct difference test has p=0.0049: these are different questions, and the latter still relies on the assumptions above. A significant result in one cohort and a nonsignificant result in another is not by itself a heterogeneity test.

The strong difference for rs28528997 is a reason to investigate its exact allele representation and cohort definition before interpreting its aligned cross-trait signals. Its rsID is multiallelic in the broader Ensembl reference. UKB and TwinGene differ in size, age distribution and FI items; separately normalizing their FIs does not make the measures identical. The UKB component includes sleep and tiredness items; removing an item was **not** tested here. Cohort relatedness/overlap has not been audited. If estimates are correlated, the standard errors above would change. No adjustment was made for winner's curse, selection, multiple candidate or analysis searches, or technical variation.

## Coding check

Each of the six variants had one exact source match by rsID, GRCh37 coordinate and unordered allele strings. Applying the inferred sign reversal to the TwinGene column and inverse-variance combining the two cohort results reproduces each existing FI meta-analysis beta within 0.0001. An expanded 97,901-variant audit found strong support for this same adapter, but it conflicts with the archive README and is **not author-confirmed**; see `twingene_direction_genomewide_audit.md`. A meta-analysis contains both components and is not an independent third cohort. Five of the six candidates have multiallelic reference rsIDs and need further study-specific validation.

The source is the [Atkins cohort summary-statistics archive](https://figshare.com/articles/dataset/Genome-Wide_Association_Study_of_the_Frailty_Index_-_Atkins_et_al_2021/9204998). The compressed source file SHA-256 is `c7ddc6e463a212c4d10ccd938970e83c9f18ecf88b4882c56f686f272a757895`; all 16,557,800 records were scanned. Full source values, reconstructed betas and comparisons are in `six_leads_fi_cohort_comparison.csv` (SHA-256 `fea96875536c1bf232555fda51344744350e1708afb18812bd0e2c2f88d02ecd`), with counts in the companion JSON. These observations do not validate a risk score or establish shared causal variants.

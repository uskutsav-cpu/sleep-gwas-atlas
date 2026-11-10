# Three exploratory clock-panel frailty leads in FinnGen insomnia

We checked the three postselected clock-panel variants highlighted in [the FI/Fried screen](stage6_hunt73_fi_ffs_screen.md) against the full FinnGen R13 insomnia GWAS (21,327,010 rows scanned; source SHA-256 `353129dc8252461a2a97087ce6de39548d0b0bf1e488a6ef36da75e9f4690049`). Using unique GRCh38 positions and allele strings returned by Ensembl, each had a matching source rsID and reference/alternate pair.

| Variant | FinnGen alternate | Insomnia OR per alternate allele | Approximate 95% CI | p |
|---|---|---:|---:|---:|
| rs2269457 | C | 1.0046 | 0.9903–1.0190 | 0.532 |
| rs2585408 | T | 1.0046 | 0.9915–1.0179 | 0.491 |
| rs10462028 | A | 1.0025 | 0.9891–1.0161 | 0.716 |

ORs and Wald intervals are calculated as exp(beta ± 1.96 × SE) from the source's alternate-allele logistic beta. All intervals cross one; this provides **no clear support for insomnia association in this dataset**, not proof of zero effects. The three were selected because of earlier frailty-only results and are not the project score variants. HUNT symptom-subtype contrasts and this clinical/register-based insomnia phenotype differ; UK Biobank/TwinGene FI, UK Biobank Fried frailty, HUNT and FinnGen also differ in sampling and ancestry. Study/relative overlap has not been fully audited. FinnGen insomnia and the project's FinnGen hospital-frailty analyses may share participants, so those two are not independent. This targeted lookup does not compare signed effects across traits, establish a shared causal variant, or validate an individual-level score.

Source: [FinnGen public results](https://www.finngen.fi/en/access_results); [row-level extraction](stage6_hunt3_finngen_insomnia_screen.csv); [streaming code](../work/hunt3_finngen_insomnia_screen.py). Ensembl GRCh38 API responses were saved in `work/hunt_rs2269457_ensembl38.json`, `work/hunt_rs2585408_ensembl38.json`, and `work/hunt_rs10462028_ensembl38.json`. No participant data were used.

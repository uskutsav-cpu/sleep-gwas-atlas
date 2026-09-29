# Draft request to the Dashti/Saxena data custodians — not sent

**Suggested contacts:** [Hassan S. Dashti, Mass General Research profile](https://researchers.mgh.harvard.edu/profile/14405948/Hassan-Dashti) (`hassan.dashti@mgh.harvard.edu` listed there); [Richa Saxena, Mass General Research profile](https://researchers.mgh.harvard.edu/profile/2993008/Richa-Saxena). Confirm the best custodian from these current institutional pages before any external contact. This file is a draft; no message or access request was sent.

**Subject:** Dashti 2019 ≥9 h versus 7–8 h UK Biobank long-sleep GWAS: variant N and released-model clarification

Dear Dr. Dashti and Dr. Saxena,

We are using the exact published long-sleep contrast (≥9 h versus 7–8 h; 34,184 cases and 305,742 controls) for a local genetic-correlation analysis. We verified that the Broad `longsumstats.txt.zip` and GWAS Catalog GCST007560 archives are byte-identical. Their ten-column file has BETA, SE, and P but no per-variant N/NMISS. We would be grateful for the smallest available clarification or data needed to represent these released statistics accurately:

1. Were `BETA_LONGSLEEP`, `SE_LONGSLEEP`, and `P_LONGSLEEP` exported directly from the primary BOLT-LMM analysis? Please provide the exact model/version, binary-phenotype treatment, any transformation of effect or test statistics, and whether a PLINK logistic sensitivity file was ever released separately.
2. Is analyzed N constant at 339,926 for **every** released variant? If not, can you share per-variant `N`/`NMISS` (or the source BOLT output field/file) matched by rsID and GRCh37 position? Per-variant cases and controls, if available, would resolve this most directly.
3. Please confirm phenotype exclusions, study-level and variant-level case/control counts, missingness handling, covariates, imputation/reference panel, variant QC, GRCh37 coordinate and effect/reference allele definitions for the released file.
4. If the original pre-portal genome-wide summary statistics contain richer fields, could those be shared with `BETA`, `SE`, `P` or signed Z and variant-level N? A checksum or release manifest would help us bind them to the public file.
5. Is a modern rerun of this **same ≥9 h versus 7–8 h** phenotype available, perhaps with REGENIE/SAIGE or another prespecified binary-trait method and per-variant N/case/control fields? A ≥10 h or 6–8 h-control analysis would be scientifically useful for a separate exploratory comparison but would not answer this exact-source question.

Our purpose is to avoid inventing a constant effective or total sample count in a downstream local analysis. A brief model/N clarification would be valuable even if no new summary file can be shared. We can follow the applicable data-use and attribution terms.

Sincerely,

Brain6 collaboration (sender/contact to be supplied by the authorized investigator)

**If only one item can be supplied:** a variant-keyed `N`/`NMISS` column for the existing file plus a written statement that the released `BETA/SE/P` are the primary BOLT-LMM output, with the exact case/control and model interpretation. That is the smallest route to a new prospective LAVA input decision; it does not predetermine that LAVA accepts the BOLT binary statistics.

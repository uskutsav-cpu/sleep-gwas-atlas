# Exact Dashti long-sleep release: source audit

The frozen [phenotype contract](phenotype_contract_v1.md) is ≥9 h habitual self-reported sleep versus 7–8 h in European-ancestry UK Biobank (34,184 cases, 305,742 controls). The study's primary association procedure was BOLT-LMM; the [published methods](https://www.nature.com/articles/s41467-019-08917-4) distinguish a separate PLINK logistic sensitivity. A coauthored [2019 sleep GWAS paper](https://www.nature.com/articles/s41467-019-11456-7) explicitly describes publicly available **long-sleep genome-wide statistics using BOLT-LMM**. This links the public analysis class more directly than the original file README, but no located export manifest states the per-variant analyzed sample counts or an exact BOLT command/transform for `longsumstats.txt`.

## Byte-level source

The existing Broad/Knowledge Portal archive and independently downloaded [GWAS Catalog GCST007560 archive](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007560/) are **byte-identical**: SHA-256 `0afd9a3ccf459d283f1b6eb6a8b8b0780f54d9d634a05f816be31171a42ec885`, 384,294,157 compressed bytes. EBI lists the file and README with last-modified timestamp 2019-05-03 12:48; that is a file timestamp, not independently verified original release date. The official README SHA-256 is `ec79d6599ad8ec0353c2d8476bd182e512b5a333519c5b6d9263bf2a9e9df0aa`.

The single ZIP member `longsumstats.txt` has SHA-256 `09e8dea3fb192f3021894e470a3130370ef152c9e5174f18e50abb2aafc36db8`, 1,094,455,329 bytes, and **14,661,601 data rows**. Exact header:

```text
SNP CHR BP ALLELE1 ALLELE0 A1FREQ INFO BETA_LONGSLEEP SE_LONGSLEEP P_LONGSLEEP
```

The actual delimiter is tab. Every row has ten fields; the full stream has zero missing-token values in each field under the audited `empty/NA/NAN/NULL/NONE/.` definition. There are 14,633,893 rsIDs and 27,708 other IDs across autosomes 1–22. The [machine-readable inventory](original_release_inventory_v1.json) contains the exact field and chromosome counts. The file has **no `N`, `NMISS`, `N_analyzed`, per-variant case/control, model, or transformation column**. A source-level missingness maximum of 10% in the paper does not prove constant per-SNP N.

## Historical Brain6 N and LAVA consequence

The historical [harmonizer](../../../scripts/01_harmonize.py) substituted a manifest-derived constant `4 × 34,184 × 305,742 / 339,926 = 122,985.40891841166` in the absence of native N. The historical QC receipt says `CONSTANT_N_EFF_FROM_MANIFEST`. The original release provided no N field to replace. [LAVA's input semantics](LAVA_SAMPLE_SIZE_SEMANTICS.md) require sample count per summary-statistic variant; input-info cases/controls set a binary proportion but do not fill or correct summary-statistic N. The historical N_eff plus observed case fraction therefore models an artificial scaled cohort. The prospective 88-locus total-N linear path [did not advance](../lava_confirmatory_pilots_v1/longsleep_linear_88_terminal_20260927.md): 44/88 TESTED, equal to canonical, with 44 NOT_RUN and no technical failures. The binary sensitivity had 47/88 TESTED but its BOLT-to-logistic assumptions remain unsupported. Neither is source-admitted.

## Search classification

The [source ledger](LONG_SLEEP_SOURCE_LEDGER.tsv) records the checked exact mirror, alternative repositories, phenotype-adjacent GWAS, and PRS-only resources. The official EBI file is a byte copy, not richer sample-size metadata. An older [Jones 2016 UKB GWAS](https://pmc.ncbi.nlm.nih.gov/articles/PMC4975467/) used the same ≥9 h versus 7–8 h contrast in only 10,102 cases and 81,204 controls; HTTP-range inspection of its [EBI GCST006685 source file](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006685/) showed BETA/SE/P but no native N column. It is an overlapping smaller precursor, not a power or N rescue. The [2023 Austin-Zimmerman UKB+MVP analysis](https://www.nature.com/articles/s41467-023-41249-y) uses **≥10 h versus 7–8 h** and is only exploratory. [SleepChart](https://www.nature.com/articles/s41586-026-10524-5) uses **>8 h versus 6–8 h**; its portal directs the genome-wide file to Synapse `syn70524096`, where an anonymous file-handle request did not permit download. [PGS019923](https://www.pgscatalog.org/score/PGS019923/) exposes PRS-CS posterior `effect_weight` values, with no variant P or SE, so it cannot substitute for GWAS summary statistics. No located source establishes an exact, richer ≥9 h versus 7–8 h rerun with per-variant N.

**Admission result:** exact phenotype available, source-valid LAVA sample-size/model representation unresolved. No new full-family run or protected promotion is justified. The [external request](LONG_SLEEP_DATA_REQUEST.md) specifies the smallest information needed to reopen this decision.

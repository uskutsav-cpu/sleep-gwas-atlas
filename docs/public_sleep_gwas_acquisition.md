# Public sleep-GWAS acquisition

This document records **source-verified public** releases selected
for Phase 0. It is an acquisition log, not an LDSC result. A downloaded file
does not become `HARMONIZATION_READY` until its real header, build,
effect-allele convention, sample-size treatment, and Phase 0 filters have been
checked.

The [Knowledge Portal Network's UK Biobank self-report page](https://www.kp4cd.org/node/235)
links directly to the five archives registered in
`config/public_gwas_sources.tsv`. The page identifies the associated papers
and states that the participants are of European ancestry. The archive member
used for napping is specified explicitly because the archive includes a
README and unrelated workspace history in addition to the summary statistics.

| Source ID | Registry traits | Publication | Public file |
| --- | --- | --- | --- |
| `dashti_2019_sleep_duration` | `sleepdur` | Dashti et al. 2019, PMID 30846698 | `sleepdurationsumstats.txt.zip` |
| `dashti_2019_short_sleep` | `shortsleep` | Dashti et al. 2019, PMID 30846698 | `shortsumstats.txt.zip` |
| `dashti_2019_long_sleep` | `longsleep` | Dashti et al. 2019, PMID 30846698 | `longsumstats.txt.zip` |
| `wang_2019_daytime_sleepiness` | `sleepiness` | Wang et al. 2019, PMID 31409809 | `Saxena.fullUKBB.DaytimeSleepiness.sumstats.zip` |
| `dashti_2021_daytime_napping` | `napping` | Dashti et al. 2021, PMID 33568662 | `Saxena_fullUKBB_Daytimenapping_summary_stats.zip` |
| `jones_2019_accelerometer_sleep` | `sleep_efficiency`, `accel_sleep_duration`, `sleep_timing` | Jones et al. 2019, PMID 30952852 | `accel_GWAS_all_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt.zip.gz` |
| `campos_2020_snoring` | `snoring` | Campos et al. 2020, PMID 32060260 | `Campos_prePMID_Snoring-mainAnalysis.gz` |
| `finngen_r9_sleep_apnoea` | `sleep_apnea` | FinnGen R9 / Kurki et al. 2023, PMID 36653562 | `finngen_R9_G6_SLEEPAPNO.gz` |

`scripts/11_materialize_public_gwas.sh` retains each downloaded archive under
the ignored `data/raw/.archives/` directory, then streams its registered
member into `data/raw/*.txt.gz`. It prints SHA-256 values after materializing,
so the local acquisition can be entered into a lab notebook without tracking
large or redistribution-restricted files in Git.

Historical short-/long-sleep aliases are not part of atlas-v1.0. The source
registry maps each selected release to exactly one locked phenotype, preventing
duplicate files from entering the Phase 1 multiple-testing family as if they
were independent traits.

## Campos 2020 snoring release

The official [GWAS Catalog record GCST009760](https://www.ebi.ac.uk/gwas/studies/GCST009760)
reports 152,302 European-ancestry cases and 256,015 controls. The primary paper
defines the binary phenotype as yes versus no to UK Biobank field 1210: whether
a partner, close relative, or friend complains about the participant's
snoring. “Don't know” and “prefer not to answer” responses were excluded.

The 258,829,716-byte original release reproduces the Catalog-supplied MD5
`c9af2918be5945ec8d8d955fdc004e14`; its project SHA-256 is
`82eeb068648959afbbccf477745820a2118f88abdff8ff0076ebb47df21f2417`.
Its literal header is `CHR BP SNP A1 A2 FREQ BETA SE P`, and `rs3094315` at
`1:752566` confirms GRCh37. Cross-checking original rows against the Catalog's
standardized file confirms that original A1 is its `effect_allele` and original
BETA is its beta on that allele. The file is therefore `SCHEMA_VERIFIED` and the
locally materialized hard link reaches `HARMONIZATION_READY`.

The real harmonization run read 11,010,158 rows and retained 7,168,629 (65.11%)
after the locked filters. Its ignored output is 156 MiB with SHA-256
`9b38e2576fdc73c1e7fbfbf8e50fd9b9ded589511e257b0b8c07ac2260f86d88`;
`data/harmonized/snoring.qc.txt` records every filter count. Snoring still
requires HapMap3 munging before `LDSC_READY`.

## FinnGen R9 sleep-apnoea release

The official R9 manifest identifies `G6_SLEEPAPNO` as 38,998 cases and 336,659
controls. Those figures replace the previous unsupported 43,901/174,054 values.
Risteys defines cases from hospital-discharge or cause-of-death ICD-10 G47.3
and ICD-9 3472 codes and excludes `G6_NARCOCATA` and `G6_SLEEPDISOTH` from
controls.

The 763,490,490-byte public archive reproduces its Google Storage MD5
`7a8e564bbab3a6e8cc3bd365874d7d60`, has SHA-256
`eb8cfc33febc044c3f608406f5f50c9be7f6973aba5effd89a0f94d1f58c1d4b`,
and passes full gzip validation. Its 20,170,208 rows use the documented header
`#chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt
af_alt_cases af_alt_controls`; FinnGen defines `alt` as the effect allele and
`beta` as its log-odds effect. The `rs2977608` coordinate is 1:832873, matching
GRCh38 and differing from GRCh37 1:768253. Source and schema are therefore
verified, but the trait remains blocked from harmonization until an explicit,
audited hg38-to-hg19 liftover is implemented.

## Jones 2019 accelerometer release

The [GWAS Catalog record GCST007803](https://www.ebi.ac.uk/gwas/studies/GCST007803)
and its release README identify a public, European-ancestry UK Biobank source
for eight accelerometer phenotypes. The archive is intentionally registered as
a `GZIP_WRAPPED_ZIP_COLUMNS` source: it is a gzip-wrapped ZIP with one
71-column member, not a misleadingly named plain-text gzip. The registered
SHA-256 and byte count refer to the downloaded outer archive.

`scripts/12_materialize_accelerometer_sleep.py` validates both containers,
requires the exact member and column layout, and streams only the three
pre-specified rank-normalized (`RAW_SIN`) phenotype groups into standardized
raw files. It does not substitute L5 or M10 timing for the actual
SPT-window-midpoint phenotype. The materialized files contain `rs6680723` at
`1:534192`; this matches the
[Ensembl GRCh37 mapping](https://grch37.rest.ensembl.org/variation/human/rs6680723?content-type=application/json),
whereas Ensembl maps GRCh37 `1:534192` to GRCh38 `1:598812`. This provides an
independent, discriminating hg19/GRCh37 sentinel, so the three continuous
accelerometer traits now pass the source/build Phase 0 curation gate.

The source pages do not replace header-level build verification. New entries
must begin with `source_status=SOURCE_PENDING` in
`config/analysis_panel.tsv` until source-registry evidence passes Phase 0. For the downloaded
Dashti 2019 sleep-duration and short-sleep files, `rs3094315` has raw
coordinate `chr1:752566`. The [Ensembl GRCh37 variation endpoint](https://grch37.rest.ensembl.org/variation/human/rs3094315?content-type=application/json)
returns `1:752566`, whereas the [GRCh38 endpoint](https://rest.ensembl.org/variation/human/rs3094315?content-type=application/json)
returns `1:817186`; this discriminating sentinel validates hg19/GRCh37 for
those two source files. The same sentinel validates the selected napping
member. The same sentinel also validates the daytime-sleepiness and long-sleep
files.
These registered files therefore reach `SOURCE_VERIFIED`. The selected
short-sleep and long-sleep rows remain blocked at `LIABILITY_H2_READY` until
their phenotype-matched population-prevalence citations are resolved; that
does not erase their source verification.

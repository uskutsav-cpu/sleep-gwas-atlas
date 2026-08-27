# Public sleep-GWAS acquisition

This document records **source-verified public** releases selected
for Phase 0. It is an acquisition log, not an LDSC result. A downloaded file
does not become `HARMONIZATION_READY` until its real header, build,
effect-allele convention, sample-size treatment, and Phase 0 filters have been
checked.

The [Knowledge Portal Network's UK Biobank self-report page](https://www.kp4cd.org/node/235)
links directly to five archives registered in `config/public_gwas_sources.tsv`.
The page identifies the associated papers and states that the participants are
of European ancestry. The UK Biobank-only insomnia and binary morning-person
releases come from the authors' CNCR share and GWAS Catalog, respectively. The
archive member used for napping is specified explicitly because the archive
includes a README and unrelated workspace history in addition to the summary
statistics.

| Source ID | Registry traits | Publication | Public file |
| --- | --- | --- | --- |
| `jansen_2019_insomnia_ukb` | `insomnia` | Jansen et al. 2019, PMID 30804565 | `Insomnia_sumstats_Jansenetal.txt.gz` |
| `jones_2019_morning_person_ukb` | `chronotype` | Jones et al. 2019, PMID 30696823 | `morning_person_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3_logORs.txt.gz` |
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

## Header and effect-column verification

`scripts/17_inspect_remote_zip_header.py` retrieves the classic ZIP central
directory and only the compressed prefix of a named member. Every request must
return HTTP 206 with the exact registered total byte count and requested
`Content-Range`; this makes header inspection reproducible without treating a
partial archive as a downloaded source artifact.

The three Dashti 2019 Catalog READMEs define `ALLELE1` as effect allele,
`ALLELE0` as reference, `A1FREQ` as effect-allele frequency, and their
trait-specific beta/SE/P columns. The selected napping archive contains its own
README defining `A1` as effect, `A2` as reference, `EAF` as A1 frequency, and
`N` as sample size. Exact registered-member headers agree with those
definitions. The binary short- and long-sleep releases are primary BOLT-LMM
betas, not the separate logistic-regression sensitivity estimates described in
the paper.

The daytime-sleepiness ZIP has no README. Its exact one-member header is
`SNP CHR BP ALLELE1 ALLELE0 A1FREQ INFO BETA SE P`. As an independent sign
check, the archive's `rs2787120` row is A/G with beta `0.00778398` and SE
`0.00137756`; those values exactly match the paper's Supplementary Data 2,
which states that effects were aligned with increasing excessive daytime
sleepiness. This verifies A (`ALLELE1`) as the beta allele for the selected
release.

The Nielsen 2018 atrial-fibrillation archive supplies a README defining
`Allele1` as effect, `Allele2` as non-effect, `Effect` as beta, and `pos` as
hg19/build37. The paper combines logistic- and Cox-model estimates and reports
relative risk as `exp(effect)`, so the schema records a generic signed log-scale
effect rather than incorrectly labeling every row as a log odds ratio.

## Nielsen 2018 atrial-fibrillation release

The exact 190,312,717-byte `AF_HRC_GWAS_ALLv11.zip` archive has SHA-256
`3ef0f55a29ba0c065df14e55cfb4fb4351ac4da90e0da28e3438426e63e4879f`.
Streaming its registered `AF_HRC_GWAS_ALLv11.txt` member into the canonical
gzip produced SHA-256
`35316d406c3188194b4efbd90de7a7bf8017fe95236f3ed7e302a8540d713649`.
The real harmonization run read 12,149,979 variants and retained 10,246,131
(84.33%). The source provides neither INFO nor allele frequency, so the QC
ledger explicitly records those unavailable filters. The harmonized contract
keeps an `FRQ=NA` placeholder; `scripts/02_munge.sh` instructs LDSC to ignore
that placeholder only when the same QC ledger proves the source field was
absent. The final HapMap3 input has 1,211,889 nonmissing signed effects after
541 allele mismatches.

The liability conversion uses `K=0.03`. The
[2016 ESC/EACTS guideline](https://academic.oup.com/eurheartj/article/37/38/2893/2334964)
estimates approximately 3% AF prevalence among adults aged 20 years or older;
the manifest labels this as an adult-population approximation rather than a
cohort-specific prevalence. The isolated five-trait checkpoint gives
liability-scale h² 0.316 (SE 0.035, Z 9.03) with intercept 1.0123, passing the
predefined gate.

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
`data/harmonized/snoring.qc.txt` records every filter count. HapMap3 munging
retained 1,179,075 nonmissing signed effects, and the liability-scale h² gate
passed in the isolated real-data checkpoint.

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
verified. The registered UCSC-chain plan now maps hg38 points to hg19, rejects
losses or ambiguous mappings, and records the chain hash in the QC ledger; the
raw source still must be locally materialized before harmonization.

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
accelerometer traits now pass the source/build Phase 0 curation gate. Their
three separate rows in `gwas_schemas.tsv` lock the exact sleep-efficiency,
SPT-duration, and SPT-midpoint `RAW_SIN` effect families; the shared archive
cannot satisfy one trait with another trait's columns.

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

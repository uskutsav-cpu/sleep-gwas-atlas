# Public sleep-GWAS acquisition

This document records the first **source-verified public** releases selected
for Phase 0. It is an acquisition log, not an LDSC result. A downloaded file
remains `TODO` until its real header, build, effect-allele convention, sample
size treatment, and Phase 0 filters have been checked.

The [Knowledge Portal Network's UK Biobank self-report page](https://www.kp4cd.org/node/235)
links directly to the five archives registered in
`config/public_gwas_sources.tsv`. The page identifies the associated papers
and states that the participants are of European ancestry. The archive member
used for napping is specified explicitly because the archive includes a
README and unrelated workspace history in addition to the summary statistics.

| Source ID | Registry traits | Publication | Public file |
| --- | --- | --- | --- |
| `dashti_2019_sleep_duration` | `sleepdur` | Dashti et al. 2019, PMID 30846698 | `sleepdurationsumstats.txt.zip` |
| `dashti_2019_short_sleep` | `shortsleep`, `shortsleep_dashti` | Dashti et al. 2019, PMID 30846698 | `shortsumstats.txt.zip` |
| `dashti_2019_long_sleep` | `longsleep`, `longsleep_dashti` | Dashti et al. 2019, PMID 30846698 | `longsumstats.txt.zip` |
| `wang_2019_daytime_sleepiness` | `sleepiness` | Wang et al. 2019, PMID 31409809 | `Saxena.fullUKBB.DaytimeSleepiness.sumstats.zip` |
| `dashti_2021_daytime_napping` | `napping` | Dashti et al. 2021, PMID 33568662 | `Saxena_fullUKBB_Daytimenapping_summary_stats.zip` |

`scripts/11_materialize_public_gwas.sh` retains each downloaded archive under
the ignored `data/raw/.archives/` directory, then streams its registered
member into `data/raw/*.txt.gz`. It prints SHA-256 values after materializing,
so the local acquisition can be entered into a lab notebook without tracking
large or redistribution-restricted files in Git.

`shortsleep`/`shortsleep_dashti` and `longsleep`/`longsleep_dashti` are
intentional aliases of the same files. The materializer uses hard links to
avoid storing duplicate data. They must be deduplicated before any Phase 1
multiple-testing family is defined; treating aliases as independent traits
would be a statistical error.

The source pages do not replace header-level build verification. New entries
must begin as `UNVERIFIED_HEADER_REQUIRED` and stay `TODO` in
`config/traits.tsv` until their actual files pass Phase 0. For the downloaded
Dashti 2019 sleep-duration and short-sleep files, `rs3094315` has raw
coordinate `chr1:752566`. The [Ensembl GRCh37 variation endpoint](https://grch37.rest.ensembl.org/variation/human/rs3094315?content-type=application/json)
returns `1:752566`, whereas the [GRCh38 endpoint](https://rest.ensembl.org/variation/human/rs3094315?content-type=application/json)
returns `1:817186`; this discriminating sentinel validates hg19/GRCh37 for
those two source files. The same sentinel validates the selected napping
member. The same sentinel also validates the daytime-sleepiness and long-sleep
files.
`sleepdur`, `sleepiness`, and `napping` are therefore the currently `CURATED`
registry rows. The short-sleep aliases remain `TODO` because their
population-prevalence citations are still missing.

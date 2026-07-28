# Public aging-GWAS acquisition

This log records source-verified public releases that have completed the
source/build portion of Phase 0. It does not make a result reproducible by
itself: the archived file still has to be materialized, harmonized, HapMap3
matched, and passed through LDSC QC locally.

| Source ID | Trait | Publication | Public file |
| --- | --- | --- | --- |
| `zenin_2019_healthspan` | `healthspan` | Zenin et al. 2019, PMID 30729179 | [`GCST007406_buildGRCh37.tsv`](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007406/GCST007406_buildGRCh37.tsv) |
| `timmers_2019_parental_lifespan` | `parental_lifespan` | Timmers et al. 2019, PMID 30642433 | [`lifegen_phase2_bothpl_alldr_2017_09_18.tsv.gz`](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009890/lifegen_phase2_bothpl_alldr_2017_09_18.tsv.gz) |
| `atkins_2021_frailty_index` | `frailty` | Atkins et al. 2021, PMID 34431594 | [`GCST90020053_buildGRCh37.tsv`](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90020001-GCST90021000/GCST90020053/GCST90020053_buildGRCh37.tsv) |

The official GWAS Catalog metadata specifies 300,447 genetically Caucasian
British UK Biobank participants and names the genome assembly as GRCh37. The
direct TSV has the documented `variant_id`, chromosome, position, effect and
other allele, effect-allele frequency, beta, standard error, and p-value
columns. Its supplied MD5 (`08f9d5861b99d0b0c00d87822eaf995c`) and the local
archive SHA-256 are recorded in `config/public_gwas_sources.tsv`.

`rs7899632` appears in the raw file as `10:100000625`, agreeing with the
[Ensembl GRCh37 variation record](https://grch37.rest.ensembl.org/variation/human/rs7899632?content-type=application/json);
the GRCh38 record is `10:98240868`. That discriminating sentinel validates
the source-build decision without a liftover.

For parental lifespan, the source README defines `a1` as effect allele, `a0`
as reference allele, and `beta1` as the log-hazard protection ratio. Its
per-SNP `n` is retained rather than silently replacing it with the registry’s
study-level count. The supplied MD5 (`9209bf29b5cf281bba35934dcf708538`) and
local SHA-256 are in the registry. Raw `rs113345124` is at `8:145793211`,
agreeing with the [Ensembl GRCh37 record](https://grch37.rest.ensembl.org/variation/human/rs113345124?content-type=application/json);
the GRCh38 position is `8:144567827`.

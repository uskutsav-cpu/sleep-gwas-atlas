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

For frailty index, the public GRCh37 TSV records `variant_id`, chromosome,
base-pair position, effect and other alleles, effect-allele frequency, beta,
standard error and p-value. The Catalog-supplied MD5
(`b8956f360155cea10d8f42b31826762e`) and local archive SHA-256 are recorded
in the registry. Raw `rs10875231` is at `1:100000012`, agreeing with the
[Ensembl GRCh37 record](https://grch37.rest.ensembl.org/variation/human/rs10875231?content-type=application/json);
the GRCh38 position is `1:99534456`.

## Verified but not eligible for the EUR panel

The public [Codd et al. Figshare release](https://figshare.com/articles/dataset/UKB_telomere_gwas_summarystats_tsv_gz/14786055)
for leukocyte telomere length was downloaded successfully (459,840,128 bytes;
supplied MD5 `f38d5d40296c8e3c8ba72357c89dc2d5`; local SHA-256
`95d5b3845c9805fcf4000ced4971040e05be5a200c77e31388c64ef33830c5a8`). Its
published header documents variant identifiers, effect and other alleles,
effect-allele frequency, beta, standard error and p-value. Raw `rs10875231`
at `1:100000012` is consistent with GRCh37.

However, the source paper's 472,174-person full-UKB analysis includes
non-European participants. The archive is registered for provenance and can
be materialized reproducibly, but `telomere_length` remains `TODO` and must
not be analysed with the EUR LDSC reference until a source-specific EUR subset
is independently identified and checked.

## Verified EUR source awaiting a citation-policy decision

The public Neale Lab round-2 left-hand grip-strength release is a GRCh37,
inverse-rank-normalized field-46 GWAS in 359,704 phenotype-complete
white-British participants. Its results archive and the associated variant
annotation passed their official manifest MD5s, and
`scripts/13_materialize_neale_grip.py` streams the two releases in exact
variant order to retain documented rsIDs, alternate effect alleles, alternate
allele frequency, beta, standard error, P, per-variant N, and INFO. The local
materialization retained 12,323,863 rsID SNVs after excluding non-rsID,
non-SNP, and source low-confidence rows. The trait-specific schema row records
the lockstep join explicitly: beta/SE/P/N come from the results file, while
rsID/position/ref/alt/frequency/INFO come from the matched annotation row.

This release is publicly citable by stable URL but was not accompanied by a
peer-reviewed article or PMID. The atlas therefore verifies its source and
schema reproducibly while keeping it blocked from `HARMONIZATION_READY` until
the project adopts an explicit citation policy for public, unpublished data
releases.

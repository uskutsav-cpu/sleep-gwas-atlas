# Public European-source substitution plan

Six locked traits cannot enter the EUR-only analysis from their currently
selected releases: three require unavailable or user-mediated access and three
are not strictly European-ancestry releases. This plan records reviewed public
substitutes without yet declaring them source-verified. Trait IDs and their
order remain unchanged.

The machine-readable candidates are in
`config/public_gwas_substitution_candidates.tsv`. Their six compressed files
total exactly 4,437,046,355 bytes (4.13 GiB). Acquisition must not begin unless
there is enough free space for the files plus harmonized and munged derivatives.

| Locked trait | Reviewed substitute | Sample | Reason and limitation |
| --- | --- | ---: | --- |
| `ms` | FinnGen R9 `G6_MS` | 2,182 cases / 373,987 controls | Public Finnish/European genome-wide data replace the unavailable IMSGC discovery meta-analysis. The substantially smaller case count and register ascertainment reduce power and must be reported. |
| `asthma` | FinnGen R9 `J10_ASTHMA_EXMORE` | 42,163 / 202,399 | Strict Finnish/European analysis replaces the public Han main-UKB file that includes all ancestries. FinnGen applies additional control exclusions, so the phenotype is not described as identical. |
| `t2d` | FinnGen R9 `T2D` | 57,698 / 308,252 | Public combined-definition T2D replaces the larger DIAMANTE file whose official route requires affirmative terms acceptance. |
| `cad` | FinnGen R9 `I9_IHD` | 63,744 / 313,533 | Strict Finnish/European data replace the predominantly-European Aragam source. The endpoint is wide ischaemic heart disease, broader than CAD, and is labeled as a proxy analysis throughout. |
| `melanoma` | FinnGen R9 `C3_MELANOMA_SKIN_EXALLC` | 2,993 / 287,137 | Public Finnish register data with all-cancer-excluded controls replace the controlled dbGaP analysis. The much smaller case count and changed control definition are explicit limitations. |
| `telomere_length` | Burren 2024 `GCST90435144` | 438,351 NFE participants | The non-Finnish-European UK Biobank stratum replaces the mixed-ancestry Codd full-UKB release. |

## Evidence and promotion gate

FinnGen's official R9 manifest supplies endpoint names, case/control counts,
and direct Google Cloud paths. The R9 summary-statistics schema is GRCh38 with
`#chrom`, `pos`, `ref`, `alt`, `rsids`, `pval`, `beta`, `sebeta`, and
`af_alt`; the registered UCSC hg38-to-hg19 chain will be required. FinnGen's
primary resource paper is Kurki et al. 2023 (PMID 36653562, DOI
10.1038/s41586-022-05473-8).

The Burren NFE file is GWAS Catalog accession `GCST90435144`, GRCh38, with
explicit chromosome, position, alleles, beta, standard error, frequency,
p-value, rsID, INFO, and per-variant N. The primary paper is Burren et al. 2024
(PMID 39192095, DOI 10.1038/s41588-024-01884-7).

A candidate can replace the selected production row only after all of these
checks pass:

1. exact byte count and upstream MD5 where supplied;
2. gzip integrity and complete-file SHA-256 registration;
3. literal header and complete row/schema audit;
4. pinned GRCh38-to-GRCh37 liftover with an explicit loss ledger;
5. phenotype-matched sample metadata and, for binary traits, a cited liability
   prevalence convention;
6. harmonization, HapMap3 munging, and predefined LDSC h2 QC without threshold
   relaxation.

Until promotion, the production manifest retains the existing selected source
and its blocker. This prevents a reviewed candidate from being mistaken for a
completed substitution.

# Full-resolution harmonization handoff

MiXeR, pleiotropic-locus discovery, fine-mapping, and TWAS require dense
post-QC GWAS rather than the HapMap3-restricted files used by LDSC. The
result-free audit identifies exactly 16 affected traits. Nine used a
memory-bounded HapMap3 source prefilter; seven used the Phase 1 HapMap3-only
rsID/coordinate map.

The seven identity-mapping traits are `parkinson`, `mdd`, `ibd`, `crohn`, `uc`,
`stroke`, and `longevity`. Their dense path is now frozen against the official
Pan-UK Biobank full variant manifest. The exact S3 object is 2,701,503,051
bytes, version `aHImXX38fFenvR2s2WKBloyhPWfsGkk7`, with single-part ETag/MD5
`e70ebc8289f762dd8d5086f54e766654`. Its exact header is independently hashed
in `config/dense_variant_map_policy.json`; the download URL also carries that
S3 version ID rather than relying on the bucket's latest object.

This resource is used only as a GRCh37 variant-identity dictionary. The builder
retains chromosome, position, reference allele, alternate allele, and rsID;
it does not transfer Pan-UKBB association results, frequency, INFO, LD, or QC
labels into another GWAS. It excludes non-autosomal/non-SNP records and every
rsID assigned to multiple distinct records. Coordinate collisions remain
explicit so lookup rejects ambiguous coordinate/allele identities rather than
selecting one.

The source is larger than the current disk can hold safely, so acquisition is
explicit and resumable:

```bash
python3 scripts/100_build_dense_variant_map.py --report-only
python3 scripts/100_build_dense_variant_map.py \
  --download --acknowledge-large-download
python3 scripts/100_build_dense_variant_map.py --build
python3 scripts/100_build_dense_variant_map.py --validate-only
```

The download requires at least 10 GiB free. The deterministic build uses a
temporary on-disk SQLite identity index, requires at least 15 GiB free, writes
fixed-mtime gzip, and publishes a SHA-256 provenance seal without overwriting
an existing map.

The complete 16-trait dispatcher derives its routes from the checksum-bound
panel, source, schema, dense-map, and liftover registries:

- seven genome-wide GRCh37 identity-map harmonizations;
- six full hg38-to-hg19 point-liftover harmonizations (`ms`, `asthma`, `t2d`,
  `cad`, `telomere_length`, and `melanoma`);
- three direct-GRCh37 GLGC harmonizations (`ldl`, `hdl`, and
  `triglycerides`).

Those three GLGC sources total exactly 6,844,892,917 bytes. On the production
host, acquire and hard-link the checksum-verified direct-gzip sources with:

```bash
bash scripts/11_materialize_public_gwas.sh --download graham_2021_ldl_eur
bash scripts/11_materialize_public_gwas.sh --materialize graham_2021_ldl_eur
bash scripts/11_materialize_public_gwas.sh --download graham_2021_hdl_eur
bash scripts/11_materialize_public_gwas.sh --materialize graham_2021_hdl_eur
bash scripts/11_materialize_public_gwas.sh --download graham_2021_triglycerides_eur
bash scripts/11_materialize_public_gwas.sh --materialize graham_2021_triglycerides_eur
```

Audit first, then materialize only on the production host:

```bash
python3 scripts/101_prepare_dense_harmonization.py --report-only
python3 scripts/101_prepare_dense_harmonization.py --materialize --all
```

The equivalent `snakemake --cores 1 dense_harmonization` DAG contains exactly
three guarded GLGC acquisition jobs, one guarded identity-map acquisition and
build, 16 trait harmonizations, and one terminal readiness gate. Its two large
download switches remain false in `config/workflow.yaml` until the transfers
are explicitly acknowledged.

Each trait is written into an isolated temporary directory, validated as a
true genome-wide route, and then published under
`data/harmonized_mixer_full/`. Existing data or QC ledgers are never replaced.
The terminal gate checks the exact trait, input path, build-normalization route,
and checksum-sealed mapping dependency, then publishes a non-overwriting
SHA-256 manifest for all 16 data/QC pairs. The dispatcher, harmonizer, mapping,
liftover, and validation scripts are hash-bound before materialization. The
detector continues to reject a QC ledger with a HapMap3-only map, even if its
old prefilter label is absent.

Pan-UKBB documents the downloadable variant manifest and GRCh37 fields on its
official downloads and per-phenotype pages. The resource citation is
Karczewski et al., *Nature Genetics* (2025), DOI
`10.1038/s41588-025-02335-7`. The identity map is a practical, versioned common
imputation-universe map; it is not a claim that every polymorphism in every
source cohort was ascertained.

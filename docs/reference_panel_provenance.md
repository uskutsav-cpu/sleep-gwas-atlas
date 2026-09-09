# EUR LD reference panel provenance

The Phase 0/1 pipeline uses the `eur_w_ld_chr` European LD-score reference
panel and its bundled `w_hm3.snplist` merge-allele list. It is the only LD
reference permitted by the repository's present EUR-only analysis plan.

| Field | Value |
| --- | --- |
| Resource | European LD scores from 1000 Genomes |
| Distribution | [Zenodo record 8182036](https://zenodo.org/records/8182036) |
| Archive | `eur_w_ld_chr.tar.gz` |
| Published archive MD5 | `e2f16343c4cfaa76caa7d0c03d26b489` |
| Verified archive SHA-256 | `9537f00eb0d163a935aaa2cf04b358b7cf21852279b9c7925802526f6060b069` |
| Archive description | Gzipped copy of the original Alkes-group `eur_w_ld_chr` distribution |
| Genome build | hg19 (GRCh37), as required by LDSC's 1000 Genomes Phase 3 EUR files |
| Installation date | 2026-07-28 |

`scripts/00_setup.sh` downloads this exact archive over HTTPS, checks the
published MD5 and repository-pinned SHA-256 before extraction, verifies all 22
autosomal LD-score files, and creates `ref/w_hm3.snplist` as a symlink to the
list inside the same archive.
The downloaded files are ignored by Git; this document, the source URL, and
the checksum are versioned instead.

## Derived GRCh37 variant-identity map

Setup also runs `scripts/19_build_hm3_variant_map.py` against the same pinned
archive. The builder joins the rsID and allele pair in `w_hm3.snplist` to the
GRCh37 chromosome and position in the 22 EUR LD-score files. It never queries
an unpinned dbSNP service and emits only the exact intersection:

| Field | Value |
| --- | ---: |
| Unique `w_hm3.snplist` rsIDs | 1,217,311 |
| Unique EUR LD-score rsIDs | 1,290,028 |
| Mapped intersection | 1,184,461 (97.3014%) |
| HapMap3 rsIDs without a bundled EUR coordinate | 32,850 |
| Output | `ref/hm3_grch37_variant_map.tsv.gz` |
| Output bytes | 9,948,937 |
| Output SHA-256 | `6775a7a0d3ca90dc74e472180b1d77103bc238129c4f969358f46307e5c306b4` |
| Validation date | 2026-08-27 |

The gzip header has a fixed timestamp and no embedded output filename, making
the map byte-reproducible across paths and repeated setup runs. Its companion
`.provenance.json` records input hashes, row counts, coverage, build, output
bytes, and output SHA-256. `scripts/01_harmonize.py` validates that sidecar and
recounts all rows before accepting the map.

`config/variant_mapping_plans.tsv` pins the output byte count and SHA-256 and
limits use to the eight selected sources whose verified schemas need identity
completion. Coordinate-based sources are
matched by GRCh37 chromosome, position, and unordered allele pair; the
coordinate-free MDD source is matched by rsID and unordered allele pair. Strand
complements are allowed, but missing, conflicting, or multiply matching rows
are dropped with explicit QC reasons. This map is deliberately an LDSC
HapMap3 gate, not a genome-wide dbSNP map and not a reference for fine-mapping.

## Pinned hg38-to-hg19 point liftover

The two source-verified GRCh38 releases are standardized with UCSC's official
[`hg38ToHg19.over.chain.gz`](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz).
UCSC documents chain filenames as `<source>To<target>` and publishes the
official MD5 in the same download directory.

| Field | Value |
| --- | ---: |
| Registered source build | GRCh38/hg38 |
| Registered target build | GRCh37/hg19 |
| Archive bytes | 1,246,411 |
| UCSC MD5 | `ff3031d93792f4cbb86af44055efd903` |
| Verified SHA-256 | `14a712e8e147d9fc8e9d87d51977b46f6f8ddb93efbe5d0843d86b6205f587b1` |
| Parsed chains | 25,374 |
| Parsed alignment blocks | 185,795 |
| Validation date | 2026-08-27 |

`config/liftover_plans.tsv` limits this operation to the verified FinnGen R9
sleep-apnoea and Bellenguez 2022 Stage I Alzheimer releases. The harmonizer
converts the source's 1-based point to chain coordinates, rejects unmapped,
non-autosomal, or multiply mapped points, and reverse-complements both alleles
when the target block is on the reverse strand. It retains the signed effect
for the correspondingly re-expressed effect allele and records the exact chain
hash and every exclusion in the trait QC ledger. This is point liftover for
summary-statistic variants, not interval or indel normalization.

The setup script downloads the chain from UCSC, checks bytes, MD5, SHA-256,
and gzip integrity, and never substitutes a newer chain silently. Users remain
responsible for UCSC's linked conditions of use and acknowledgment guidance.

This provenance does not validate a GWAS. Before any LDSC run, the GWAS source
must independently be confirmed as EUR and hg19 (or have a separately
documented build decision); the panel must never be used to excuse an ancestry
or build mismatch.

## LDSC Python 3 compatibility

The reproducible installer clones CBIIT/ldsc's `ldsc39` branch and applies the
version-checked patch at `patches/ldsc39-python3-compressed-header.patch`. At
the reviewed revision, `munge_sumstats.py` opens a `.gz` or `.bz2` header in
binary mode but strips a text newline, which fails under Python 3 before any
summary statistic is parsed. The patch changes only that header stream to text
mode; it is not a statistical-method modification. Setup stops if the patch
does not match the checked-out LDSC source, preventing an unreviewed upgrade.

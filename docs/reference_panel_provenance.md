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
| Archive description | Gzipped copy of the original Alkes-group `eur_w_ld_chr` distribution |
| Genome build | hg19 (GRCh37), as required by LDSC's 1000 Genomes Phase 3 EUR files |
| Installation date | 2026-07-28 |

`scripts/00_setup.sh` downloads this exact archive over HTTPS, checks the
published MD5 before extraction, verifies all 22 autosomal LD-score files, and
creates `ref/w_hm3.snplist` as a symlink to the list inside the same archive.
The downloaded files are ignored by Git; this document, the source URL, and
the checksum are versioned instead.

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

# Phase 1 methods

## Software and reference data

**LDSC**: maintained CBIIT Python-3.9 fork (`ldsc39` branch), not the stale
Python 2.7 original. Reachable at `ldsc/` via symlink to a shared install;
`.ldsc-env/` is its Python 3.9.23 environment.

**LD reference**: `eur_w_ld_chr` (1000 Genomes Phase 3, European, HapMap3
SNPs), 48 files. `w_hm3.snplist` is the `--merge-alleles` target. EUR panel
against EUR summary statistics throughout — no ancestry substitution.

## Source selection

Every source located through the GWAS Catalog REST API by verified PMID or by
trait, never by guessing a URL. Where a publication maps to several accessions
the **European-only stratum** was selected deliberately, and the multi-ancestry
alternative recorded. Every PMID was confirmed against NCBI E-utilities before
use.

## Provenance and reproducibility

`config/public_gwas_sources.tsv` records per trait: accession, PMID, landing
page, direct URL, archive member, access status, declared vs **verified**
build, byte count, SHA-256, and variant count.

Raw and harmonized files are deleted after munging (disk is single-digit GB);
the URL and hash make each input exactly reproducible. This was validated: four
traits re-downloaded during this work returned **byte-identical SHA-256s**, and
four sleep traits independently matched hashes recorded by a separate
workstream.

## Genome build verification

Build is **proven from file contents, never from filenames**. `config/build_anchors.tsv`
holds nine SNPs whose GRCh37 and GRCh38 positions differ by ≥ 1 kb, fetched
from Ensembl on both assemblies with their allele strings.

A finding worth recording: **position alone does not discriminate build** in a
densely imputed file, because a variant exists at both builds' coordinates. The
check therefore requires position **and** allele agreement, and reports LOW
CONFIDENCE below three matched anchors.

## Harmonization

`scripts/01_harmonize.py`, filters tagged `[CDG3]` (from the Nature methods) or
`[LDSC]` (tool requirement): strand-ambiguous removed, INFO > 0.6, MAF > 1%,
MHC chr6:25–34 Mb excluded, duplicates removed, SNP-specific N ≥ 50% of median.
Every filter's drop count and the retention percentage go to a per-trait QC
ledger (`data/harmonized/<trait>.qc.txt`).

**Chunked streaming is used for large files** (`--chunksize`). On an 8 GB
machine the single-pass path is killed outright by 13 M-row inputs. Chunked and
in-memory paths were regression-tested to produce identical output (BMI:
1,973,592 SNPs kept either way). Two documented compromises in chunked mode:
duplicate rsIDs are removed within chunks (LDSC munge removes cross-file
duplicates), and the SNP-specific N filter uses the first chunk's median.

## Heritability and the QC gate

h² is estimated on the **observed scale**. Every population prevalence in the
registry is `UNKNOWN` because none carried a citation, and rule 1 holds that an
uncited prevalence is worse than a blank. Liability conversion is deferred.

Gate: **Z = h²/SE ≥ 4** (below this rg is not interpretable) **and LDSC
intercept ≤ 1.20** (above this suggests confounding rather than polygenicity).
These are different failure modes and are never reported as the same thing.

**N_eff × h² > 12,000** (CDG3's MiXeR inclusion threshold) is recorded
separately and is *not* part of the rg eligibility gate.

Because the conversion rescales h² and its SE by the same constant, **Z is
invariant to the liability choice** — confirmed empirically here (shortsleep
23.7 → 23.6, longsleep 13.6 → 13.8).

## Genetic correlation and multiple testing

Bivariate LDSC. The **primary family is declared in advance**: every eligible
core sleep/circadian trait × every eligible disease trait. BH-FDR is applied
across that full primary family, not chosen after seeing results. Secondary
sleep and aging traits are held in a separate panel so they cannot inflate the
primary family.

Cross-trait intercepts (`gcov_int`) are retained in every output row; sample
overlap is expected wherever UK Biobank appears on both sides and is handled by
that intercept rather than assumed absent.

## Known non-independence

`sbp` / `dbp` / `pulse_pressure` are one cohort of 1,028,980, and pulse
pressure is by definition SBP − DBP. Lipid and glycaemic traits likewise share
cohorts. These are curated separately but must not be counted as independent
tests.

# LAVA local-architecture workflow

## Locked design

The atlas uses LAVA 0.1.5 at official commit
`e729a245f7b6923967a96804fbf5246eadf2d6c6`. Its bundled GRCh37/hg19
partition contains 2,495 loci and is checksum locked in
`config/lava_analysis_policy.json`. The analysis includes all 45 locked traits
and all 396 sleep-by-non-sleep pairs. Global LDSC correlation magnitude or
significance never filters local tests.

The local-univariate family is fixed at 45 x 2,495 = 112,275 planned tests.
The two traits in a locus must each pass the preregistered Bonferroni threshold
(`0.05 / 112275`) before that locus-pair enters bivariate LAVA. Benjamini-Hochberg
FDR is then computed over the complete family of successfully tested eligible
locus-pair rows. Phase-1 QC-failed traits remain labelled sensitivity analyses
rather than being silently removed.

Sample overlap is estimated from the completed 45-trait cross-trait LDSC
intercept matrix, standardized exactly as the official LAVA tutorial specifies.
The resulting correlation matrix is positive definite (minimum eigenvalue
0.0270711; condition number 83.55).

## Software and inputs

Run the checksum-pinned software setup and input preparation with:

```bash
bash scripts/30_setup_lava.sh
python3 scripts/31_prepare_lava.py
```

`30_setup_lava.sh` installs exact CRAN sources for `keep` 1.0,
`matrixsampling` 2.0.0, and `cpp11` 0.5.2 before installing the pinned LAVA
commit. On Apple Silicon it uses the system Apple Clang through
`environment/macos-arm64.Makevars`. The setup extracts the official locus file
but deliberately does not download an LD reference.

`31_prepare_lava.py` verifies all munged files contain `SNP`, `A1`, `A2`, `Z`,
and `N`; checks binary case, control, and prevalence metadata; constructs the
LAVA input-info table; validates and writes the overlap matrix; emits the exact
396-pair manifest; and fingerprints every summary-statistics file.

## Reference storage blocker

LAVA 0.1.5 strongly recommends its UK Biobank LD v1.1 reference for European
analyses and warns that the older 1,000 Genomes reference can inflate type-I
error and bias local heritability. The seven official compressed archives total
14,110,596,095 bytes (13.14 GiB), and the official documentation reports 15 GiB
uncompressed. The downloader therefore requires at least 35 GiB free so the
archives, extracted files, and a safety margin can coexist.

Inspect the exact URLs, byte counts, target, and available storage without
downloading:

```bash
bash scripts/32_download_lava_reference.sh
```

The actual download is separately opt-in:

```bash
bash scripts/32_download_lava_reference.sh --download
```

As of 28 August 2026, the local volume has only roughly 4–6 GiB free, so the
production reference has not been downloaded. No 1,000 Genomes substitution is
permitted merely to fit the current disk.

## Production and validation

Once the complete reference is present, run:

```bash
snakemake --cores 1 lava
python3 scripts/34_validate_lava.py
```

The production runner checkpoints each locus under
`results/checkpoints/lava/`. It writes the canonical tables only after all
2,495 checkpoints exist with the same input/reference fingerprint. Every
locus-trait receives a row, including explicit processing or phenotype-drop
statuses, and every locally eligible locked pair is either tested or records a
failure. The validator enforces coverage, family sizes, Bonferroni values,
eligible-pair identity, FDR, version, reference, and predefined failure limits.

Official sources:

- <https://github.com/josefin-werme/LAVA>
- <https://github.com/josefin-werme/LAVA/blob/master/REFERENCE.md>

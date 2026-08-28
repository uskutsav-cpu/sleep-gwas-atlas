# MiXeR polygenic-overlap workflow

## Locked design

The atlas pins the latest non-prerelease MiXeR release, 2.2.1, at tag commit
`cf65c57d5d1ad76597db1d4fa3907f1d711d39e7`. The Linux/amd64 container is
locked to manifest
`sha256:5bf54ddd6f7f81b93eeb5450b1a6dc809d1fcf926b3a815d6d514f36ce04f51c`;
its compressed layers total 2,106,984,629 bytes (1.96 GiB). The exact source,
documentation, and reference-data commits observed when the policy was frozen
are recorded in `config/mixer_analysis_policy.json`.

Every trait is run through 20 `fit1`/`test1` replicates using the official 20
random approximately 600K-SNP extract sets. A trait enters bivariate MiXeR only
when all replicates and combined outputs exist, model parameters are finite,
and the combined univariate fit has AIC > 0. This follows the official FAQ:
negative AIC indicates insufficient information to justify the MiXeR model.
Positive AIC with non-positive BIC is retained but labelled
`BORDERLINE_POWER_AIC_ONLY`. Every sleep-by-non-sleep pair for which both traits
pass is then run through 20 `fit2`/`test2` replicates. Global LDSC correlation
does not filter MiXeR pairs.

## Input contract

MiXeR 2.2.1 accepts LDSC-like `SNP, CHR, BP, A1, A2, N, Z` summary statistics,
but its current guidance explicitly advises against HapMap3-filtering and says
older HapMap3-constrained inputs should be regenerated. Consequently, the nine
traits that used a memory-saving HapMap3 source prefilter for LDSC cannot be
promoted to primary MiXeR inputs:

```text
ms, asthma, t2d, ldl, hdl, triglycerides, cad, telomere_length, melanoma
```

Six still have their full registered raw archives locally; the three lipid
traits require their full source archive to be rematerialized. Full MiXeR
harmonizations belong under `data/harmonized_mixer_full/` with matching QC
ledgers and no `HAPMAP3_RSID_ALLOWLIST` provenance.

Audit the exact input and host status without creating large files:

```bash
python3 scripts/35_mixer_preflight.py --report-only
python3 scripts/36_prepare_mixer_inputs.py
```

Only after all 45 full inputs are available and the storage plan has been
reviewed should the conversion be requested:

```bash
python3 scripts/36_prepare_mixer_inputs.py --materialize
```

The converter streams each GWAS, calculates `Z=BETA/SE`, writes deterministic
gzip, and records input/output checksums and row counts. It fails before writing
if any trait remains HapMap3-prefiltered or if free space is below 1.5 times the
compressed sources plus 2 GiB.

## Required analysis host

Official MiXeR containers are x86-only. The real-data tutorial specifies at
least 32 GB RAM and recommends at least 16 physical cores; the official
reference repository is about 14 GB. The current laptop is Apple M1/arm64 with
8 GB RAM, 8 logical cores, and only about 5.4 GiB free. It therefore fails
architecture, memory, CPU, storage, reference, and full-input preflights.
Unsupported amd64 emulation is not treated as a scientific production run.

On a supported x86 host, materialize the official reference at
`ref/mixer/reference/`, preserving reference commit
`a4104bf34ed0509a6daa5b06594c40f0b655871c`. Inspect the 1.96 GiB image plan,
then explicitly pull it:

```bash
bash scripts/37_pull_mixer_image.sh
bash scripts/37_pull_mixer_image.sh --pull
```

The first command never downloads. Reference acquisition is also deliberately
not automated because the official repository currently warns that its Git LFS
quota is broken and routes users to an external Dropbox copy.

## Production tasks and validation

The task runner maps directly to the official 2.2.1 real-data sequence:

```bash
bash scripts/38_run_mixer_task.sh univariate insomnia 1
bash scripts/38_run_mixer_task.sh combine-univariate insomnia
bash scripts/38_run_mixer_task.sh bivariate insomnia bmi 1
bash scripts/38_run_mixer_task.sh combine-bivariate insomnia bmi
```

Replicate indices are 1–20 and seeds are fixed at 1001–1020. Production should
submit these as cluster arrays. The runner refuses to proceed unless the whole
preflight passes and the pinned image already exists, so it never triggers an
implicit multi-gigabyte pull.

After combining every univariate model, and then every pair implied by those
univariate results, publish and validate canonical tables:

```bash
python3 scripts/39_collate_mixer.py --univariate-only
python3 scripts/39_collate_mixer.py
python3 scripts/40_validate_mixer.py
```

The validator enforces 45 ordered univariate models, 20-replicate provenance,
the AIC/BIC power rules, the exact eligible sleep-by-non-sleep pair family,
finite parameters, plausible overlap/concordance ranges, and the pinned image
and reference identifiers.

Official sources:

- <https://github.com/precimed/mixer>
- <https://github.com/precimed/gsa-mixer/releases/tag/v2.2.1>
- <https://github.com/comorment/mixer>

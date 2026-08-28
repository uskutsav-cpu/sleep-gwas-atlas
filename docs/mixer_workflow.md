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
older HapMap3-constrained inputs should be regenerated. Consequently, 16
traits whose LDSC artifacts used either a memory-saving HapMap3 source
prefilter or a HapMap3-only variant-identity map cannot be promoted to primary
MiXeR inputs:

```text
ms, asthma, t2d, ldl, hdl, triglycerides, cad, telomere_length, melanoma,
crohn, ibd, uc, mdd, longevity, parkinson, stroke
```

Thirteen still have their full registered raw archives locally; the three lipid
traits require their full source archives to be rematerialized. The seven
variant-map cases now have a result-free genome-wide GRCh37 identity contract
based on the exact Pan-UKBB variant manifest. Its guarded builder and the exact
16-trait routing are documented in `docs/dense_harmonization_workflow.md`.
Full MiXeR harmonizations belong under `data/harmonized_mixer_full/` with
matching QC ledgers and neither `HAPMAP3_RSID_ALLOWLIST` nor HapMap3-only map
scope provenance.

Audit the exact input and host status without creating large files:

```bash
python3 scripts/35_mixer_preflight.py --report-only
python3 scripts/100_build_dense_variant_map.py --report-only
python3 scripts/101_prepare_dense_harmonization.py --report-only
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
compressed sources plus 2 GiB. A successful conversion also writes an immutable
companion lock binding the ordered 45 traits, panel, policy, manifest, and every
converted file identity. Existing inputs or locks are never silently replaced.

## Required analysis host

Official MiXeR containers are x86-only. The real-data tutorial specifies at
least 32 GB RAM and recommends at least 16 physical cores. The broader official
reference repository is larger, but the exact 64-file family consumed by this
analysis is 6,579,093,199 bytes: 22 BIM files, 22 run4 LD files, and 20
replicate-specific extract lists. `config/mixer_reference_files.tsv` pins every
payload to the Git LFS SHA-256 and pointer identity at official `comorment/mixer`
commit `a4104bf34ed0509a6daa5b06594c40f0b655871c`. The current laptop is
Apple M1/arm64 with 8 GB RAM, 8 logical cores, and only about 2.0 GiB free. It
therefore fails architecture, memory, CPU, storage, reference, and full-input preflights.
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
quota is broken and routes users to an external Dropbox copy. After acquisition,
hash all 6.58 GB once and seal the exact local family before any fit:

```bash
python3 scripts/35_mixer_preflight.py --seal-reference --report-only
```

Ordinary task preflights then compare the immutable seal and every live file
size without rehashing 6.58 GB for each array job. Final result validation
rehashes all returned task artifacts.

## Production tasks and validation

On the supported host, freeze the complete univariate family before execution:

```bash
python3 scripts/mixer_tasks.py univariate
python3 scripts/mixer_tasks.py univariate --write
```

The first command is a read-only plan. The second is allowed only after the
whole host, reference, source, and converted-input preflight passes. It creates
900 replicate tasks plus 45 combine tasks, with exact commands, seeds, expected
outputs, and a pre-result lock. The task runner refuses commands absent from
that lock and maps each row directly to the official 2.2.1 sequence:

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

After combining every univariate model, publish and validate it before the
bivariate family can be defined. Then freeze exactly the pairs whose two
univariate models pass AIC eligibility:

```bash
python3 scripts/39_collate_mixer.py --univariate-only
python3 scripts/40_validate_mixer.py --univariate-only
python3 scripts/mixer_tasks.py bivariate
python3 scripts/mixer_tasks.py bivariate --write
python3 scripts/39_collate_mixer.py
python3 scripts/40_validate_mixer.py
```

The validator enforces 45 ordered univariate models, 20-replicate provenance,
the AIC/BIC power rules, the exact eligible sleep-by-non-sleep pair family,
finite parameters, plausible overlap/concordance ranges, and the pinned image
and reference identifiers. A scientifically complete zero-eligible-pair family
is represented by a schema-bearing empty bivariate table instead of being
misreported as a pipeline failure. Each phase receives immutable provenance
binding its task lock, exact input/reference seals, every returned artifact,
and canonical table. This permits transfer by encrypted disk or `rsync`; no Git
push is required.

The same sequence is represented by the `mixer` Snakemake target. It converts
all 45 full inputs, freezes and runs the 900 univariate replicates plus 45
combines, checkpoints the observed eligibility family, then schedules exactly
the eligible bivariate replicates and combines before final collation:

```bash
snakemake --cores 16 mixer
```

The 64-file reference must still be staged manually and sealed because the
official repository currently routes around a broken Git LFS quota. The image
is accepted when already present; otherwise
`acknowledge_mixer_container_pull: true` must be set before Snakemake may pull
the exact 2,106,984,629-byte image.

Official sources:

- <https://github.com/precimed/mixer>
- <https://github.com/precimed/gsa-mixer/releases/tag/v2.2.1>
- <https://github.com/comorment/mixer>

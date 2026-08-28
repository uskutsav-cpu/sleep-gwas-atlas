# Pleiotropic-locus workflow

## Locked methods and scope

The primary shared-locus layer uses two complementary methods over the same
locked 396 sleep-by-non-sleep pairs: PLACO+ 0.2.0 and conjunctional FDR from
the official `precimed/pleiofdr` implementation. Genome-wide LDSC correlation
does not select pairs. T2D and melanoma remain labelled QC-failed sensitivity
pairs, but they are not silently removed from the scan family.

PLACO+ is pinned to commit
`3ba3cae1d323ad117fb4540e620bcefa79f70663` and source SHA-256
`fb684a8ed88f27dd138f5c2e8613b092904e84f364030d036058f17db95b7124`.
The nuisance variance/correlation estimates use marginal P > 1e-4 variants,
and the primary numerical path excludes variants where either squared Z score
exceeds 80. The variant threshold is 5e-8/396 across the frozen pair family.

ConjFDR is pinned to official commit
`0da963cac22fe8de9030166d7aea974bcbb0a367`, uses 500 random-pruning
iterations, removes the extended MHC while fitting, corrects sample overlap,
and uses the standard conjFDR <= 0.05 discovery threshold. A canonical shared
locus must be supported by both PLACO+ and conjFDR in the same one of the
official 2,495 approximately independent GRCh37 blocks. Method-specific rows
remain available so agreement is auditable.

Statistical pleiotropy is not interpreted as proof of a shared causal variant,
biological mechanism, or causal direction. Those claims require the later
fine-mapping and colocalization layers.

The canonical locus table retains the two signed PLACO+ Z statistics at its
lead variant and reports their effect pattern as `+/+`, `+/-`, `-/+`, or
`-/-`. This direction label is descriptive on the shared reference allele; it
is not a causal-direction estimate.

## Current preflight

The official pleioFDR reference is exactly 2,383,912,974 bytes (2.22 GiB) and
the upstream workflow requires MATLAB and at least 16 GB RAM. It is not
downloaded implicitly. This Apple M1 laptop has 8 GB RAM, no MATLAB, about
2.2 GiB free, and only 36/45 full post-QC inputs. The nine inputs that were
source-prefiltered to HapMap3 for LDSC block 108 of 396 pair scans:

```text
ms, asthma, t2d, ldl, hdl, triglycerides, cad, telomere_length, melanoma
```

Run the read-only audits with:

```bash
bash scripts/41_setup_pleiotropy.sh
python3 scripts/42_pleiotropy_preflight.py --report-only
python3 scripts/43_prepare_pleiotropy_pairs.py --report-only
python3 scripts/44_materialize_pleiotropy_pair.py insomnia__bmi
```

The fourth command is also read-only by default: it validates the locked
pair, source hashes, policy hash, and available scratch space. Add
`--materialize` only on a production host with enough free space; the output
is written deterministically and published atomically after full allele and
coordinate checks.

The pair planner is likewise non-mutating under `--report-only`. It publishes
the immutable 396-row manifest only when all 396 pairs have two full inputs;
an incomplete planning snapshot cannot authorize partial production tasks.

After materialization, freeze and run one PLACO+ scan with:

```bash
python3 scripts/45_prepare_placo_task.py insomnia__bmi
.r-env/bin/Rscript scripts/46_run_placo_pair.R \
  results/pleiotropy/tasks/insomnia__bmi.tsv \
  results/pleiotropy/tasks/insomnia__bmi.lock.tsv --execute
```

The runner validates the task, pair input, and pinned upstream source before
loading results. It evaluates the full variant family, records all numerical
failures, fails if their fraction exceeds the locked limit, and retains both
the conventional and stricter 396-scan-family indicators. A real pair scan is
substantial production compute and is never launched by a read-only preflight.
Every task is also bound to the exact task builder and R runner, and existing
task or result artifacts are never silently replaced.

ConjFDR also needs the official 9,545,380-variant ordering template. Its exact
274,423,819-byte payload is independently pinned by SHA-256. It is not
downloaded by default:

```bash
bash scripts/41_setup_pleiotropy.sh --download-template
python3 scripts/47_prepare_pleiofdr_trait.py insomnia
```

The trait command audits only unless `--materialize` is supplied. A real
materialization requires the pinned SciPy environment and at least 8 GiB free;
it writes dense `logpvec`/`zvec` vectors in exact reference order, with missing
variants represented as `NaN` as required by the upstream implementation.

Once both trait MAT files and the 2.22 GiB LD reference exist, freeze and run
one exact conjunction-FDR task with:

```bash
python3 scripts/48_prepare_conjfdr_task.py insomnia__bmi
python3 scripts/49_run_conjfdr_pair.py \
  results/pleiotropy/conjfdr_tasks/insomnia__bmi.tsv \
  results/pleiotropy/conjfdr_tasks/insomnia__bmi.lock.tsv --execute
```

The runner copies the clean pinned checkout into ignored scratch space and
applies the checksum-locked one-line option patch that enables upstream's
Mahalanobis sample-overlap correction. It refuses experimental Octave, less
than 16 GiB RAM, less than 20 GiB free disk, drifted inputs, or execution
without `--execute`.

Once the exact PLACO+ source, clean pleioFDR commit, variant template, and LD
reference are all present, setup hashes the complete family once and writes
`ref/pleiofdr/runtime.provenance.json`. Pair tasks verify that immutable seal
and live reference size instead of rereading 2.22 GiB for every array job. The
final collator rehashes the LD reference once before publication and validates
every PLACO+ task, conjunction-FDR task, returned CSV/MAT/log, and runner hash.

After all 396 pair completions exist, the only canonical publisher is:

```bash
python3 scripts/50_collate_pleiotropy.py
```

It validates every task and result hash, assigns method-specific hits to the
locked 2,495 blocks, and publishes `results/atlas/shared_loci.tsv` only for
pair/block combinations supported by both methods. `--report-only` is a
non-mutating completeness audit; `--validate-only` recomputes and byte-checks
an existing canonical result. Canonical outputs and their provenance are
immutable and include checksums for the full 396-pair returned-artifact family.

The equivalent opt-in Snakemake targets are `placo_pair`, `conjfdr_pair`, and
`pleiotropy`. They are intentionally outside the default `all` target because
they materialize large files and launch production-scale computation.

The setup script only prints the software/reference plan by default. The small
code downloads require `--download-software`. The 2.22 GiB reference requires
the separate `--download-reference` flag and a 20 GiB free-space preflight.
No download has been requested on the current laptop.

Primary sources:

- <https://github.com/RayDebashree/PLACO>
- <https://github.com/precimed/pleiofdr>

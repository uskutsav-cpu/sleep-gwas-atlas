# sleep-gwas-atlas

Phase 0–1 pipeline for the Sleep/Circadian Genetic Atlas: harmonize GWAS
summary statistics, estimate SNP heritability, apply a power QC gate, run
pairwise genetic correlations, and produce the Figure 2 heatmap.

Built to mirror the Methods of Grotzinger/Werme et al., *Nature* 649:406–415
(2026). See `methods_map.md`.

## Status

| Component | State |
|---|---|
| `01_harmonize.py` | tested on synthetic sumstats (tab + space delimited, BETA + OR formats) |
| `05_collate.py` | tested on synthetic LDSC logs; BH-FDR verified monotone |
| `06_heatmap.py` | tested via the reproducible smoke test below; provenance footer + synthetic-input guard |
| `00_setup.sh`, `02_munge.sh`, `03_h2_qc.sh`, `04_rg.sh` | **written but NOT executed** — they need conda + LDSC + real reference panels, which the machine they were authored on could not reach |

Treat the four shell scripts as a first draft that will need a debugging pass
on real files. That is normal and expected; the Python is the part that has
actually been exercised.

## Run order

```bash
bash scripts/00_setup.sh                     # once: LDSC + conda env + ref panels
# put raw sumstats in data/raw/, named as in config/traits.tsv
bash scripts/02_munge.sh insomnia mdd        # harmonize + munge
bash scripts/03_h2_qc.sh insomnia mdd        # h2 + QC gate  <- the decision point
# edit config/traits.tsv: set status=PASS for survivors
bash scripts/04_rg.sh                        # all PASS sleep x all PASS disease
```

To see the whole thing move without any real data:

```bash
python3 scripts/make_test_data.py --out data/raw
python3 scripts/01_harmonize.py --trait insomnia --config config/traits.tsv \
    --infile data/raw/insomnia.txt.gz --outdir data/harmonized
cat data/harmonized/insomnia.qc.txt          # the auditable QC ledger
```

## Smoke test (collate + gate + figure, no LDSC required)

Synthetic LDSC logs, so the QC gate and the plotting code can be exercised
without conda or reference panels. Output goes to `results/_smoketest/` and
never to `results/logs/` or `results/tables/` — fake and real LDSC output must
not be able to share a directory, because `05_collate.py` globs a directory
and cannot tell them apart.

```bash
python3 scripts/make_test_data.py --logs
python3 scripts/05_collate.py --mode h2 --logdir results/_smoketest --out results/_smoketest/h2_summary.tsv
python3 scripts/05_collate.py --mode rg --logdir results/_smoketest --out results/_smoketest/rg_matrix.tsv
python3 scripts/06_heatmap.py --rg results/_smoketest/rg_matrix.tsv --config config/traits.tsv \
    --out results/_smoketest/fig2_smoketest.png --provenance-label "SYNTHETIC — NOT REAL DATA"
```

What it deliberately covers: negative rg (the blue half of the colorbar),
cells that survive BH-FDR at every star level *and* four that do not
(blank-star path), a trait that fails the gate on **power** (`longsleep`,
Z=2.0) and a separate one that fails on **confounding** (`shortsleep`,
intercept 1.35). Values are hardcoded, never jittered.

`06_heatmap.py` refuses to render an rg table from `_smoketest/` unless
`--provenance-label` is passed, and the label becomes a diagonal watermark.
Every figure — real or synthetic — also gets a footer with its source table's
path, mtime, and content hash, so a figure can always be traced back to the
numbers that made it.

**Generated output is gitignored and is never committed.** A figure that
outlives its inputs is an orphan, and an orphan figure with plausible numbers
on it is indistinguishable from a result. This repo has already produced one.

## Things you must fix before this is real

1. **`config/traits.tsv` sample sizes and population prevalences are
   placeholders.** I filled them from memory of the literature. Every one
   needs checking against the actual paper. A wrong `pop_prev` silently
   corrupts every liability-scale h² downstream, and it is exactly the kind
   of number a mentor spot-checks.
2. **Verify the insomnia release you download excludes the 23andMe component**
   unless you have access.
3. **Confirm the LD-score download URLs resolve.** Broad's hosting has moved;
   both fallbacks are in `00_setup.sh`, and the LDSC wiki is authoritative.
4. **Decide hg19 vs hg38 with your mentor** before harmonizing anything. CDG3
   used GRCh37 throughout.

## The QC gate

`03_h2_qc.sh` → `05_collate.py --mode h2` classifies each trait PASS/DROP on:

- **Z = h²/SE ≥ 4** — below this, rg estimates are not interpretable
- **LDSC intercept ≤ 1.2** — above this suggests confounding, not polygenicity

Expect `longsleep` and possibly `shortsleep` to fail. They are dichotomised
tails of a continuous phenotype and carry little heritability. That result is
not a setback; it is the finding that scopes Phase 4, and it belongs at the
top of your next meeting agenda.

## Layout

```
config/traits.tsv        trait registry = the Phase 0 metadata table
scripts/01_harmonize.py  QC with an auditable ledger per trait
scripts/05_collate.py    LDSC logs -> tidy tables + QC gate + BH-FDR
scripts/06_heatmap.py    Figure 2
data/harmonized/*.qc.txt every SNP dropped, and why
results/tables/          h2_summary.tsv, rg_matrix.tsv
results/figures/         fig2_rg_heatmap.png
```

Each QC threshold in `01_harmonize.py` is tagged `[CDG3]` where it comes from
the Nature paper's Methods or `[LDSC]` where the tool requires it. Read those
tags before the meeting — you will be asked why you filtered at INFO 0.6, and
"it was in the script" is not the answer you want to give.

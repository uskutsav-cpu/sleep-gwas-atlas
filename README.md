# sleep-gwas-atlas

Phase 0–1 pipeline for the Sleep/Circadian Genetic Atlas: harmonize GWAS
summary statistics across a **90-trait catalog**, estimate SNP heritability ($h^2$),
apply a power QC gate, evaluate a **dual-strategy ("use both") comparison framework**,
run pairwise genetic correlations ($r_g$), and produce the Figure 2 heatmap.

Built to mirror the Methods of Grotzinger/Werme et al., *Nature* 649:406–415
(2026). See `methods_map.md`.

## Dual-Strategy Framework ("Use Both")

Rather than arbitrarily picking one definition or dataset release, the pipeline
evaluates both side-by-side:

1. **Short/Long Sleep Phenotype Definitions**:
   - Dashti 2019 (<7h / >=9h, Ncase=106k/34k) **vs** Austin-Zimmerman 2023 (<=5h / >=10h, extreme tail).
2. **Insomnia Summary Statistics**:
   - Jansen 2019 UKB-only (Public, N=386k) **vs** Jansen 2019 UKB+23andMe (Full, N=1.33M).
3. **Liability Scale & Ascertainment**:
   - Standard sample size $N_{total}$ **vs** CDG3 effective sample size $N_{eff}$ ($N_{eff} \times h^2 > 12,000$ MiXeR threshold).

`05_collate.py` automatically generates a dedicated comparative summary table:
`results/tables/problem_comparison_summary.tsv` (or `results/_smoketest/problem_comparison_summary.tsv` during smoke testing).

## Status

| Component | State |
|---|---|
| `01_harmonize.py` | tested on synthetic sumstats (tab + space delimited, BETA + OR formats, build & $N_{eff}$ logging) |
| `05_collate.py` | tested on synthetic LDSC logs; BH-FDR verified monotone; generates dual-strategy comparison table |
| `06_heatmap.py` | tested via the reproducible smoke test below; dynamic layout scaling for 90 traits; provenance footer + synthetic-input guard |
| `00_setup.sh`, `02_munge.sh`, `03_h2_qc.sh`, `04_rg.sh` | **written but NOT executed** — they need conda + LDSC + real reference panels |

## Run order

```bash
bash scripts/00_setup.sh                     # once: LDSC + conda env + ref panels
# put raw sumstats in data/raw/, named as in config/traits.tsv
bash scripts/02_munge.sh insomnia mdd        # harmonize + munge
bash scripts/03_h2_qc.sh insomnia mdd        # h2 + QC gate  <- the decision point
# edit config/traits.tsv: set status=PASS for survivors
bash scripts/04_rg.sh                        # all PASS sleep x all PASS disease
```

To see the 90-trait pipeline move without any real data:

```bash
python3 scripts/make_test_data.py --logs
python3 scripts/05_collate.py --mode h2 --logdir results/_smoketest --out results/_smoketest/h2_summary.tsv
python3 scripts/05_collate.py --mode rg --logdir results/_smoketest --out results/_smoketest/rg_matrix.tsv
python3 scripts/06_heatmap.py --rg results/_smoketest/rg_matrix.tsv --config config/traits.tsv \
    --out results/_smoketest/fig2_smoketest.png --provenance-label "SYNTHETIC — NOT REAL DATA"
```

## The QC gate & Dual-Strategy Table

`05_collate.py --mode h2` classifies each trait PASS/DROP on:

- **Z = h²/SE ≥ 4** — below this, rg estimates are not interpretable
- **LDSC intercept ≤ 1.2** — above this suggests confounding, not polygenicity

Additionally, it outputs `problem_comparison_summary.tsv` comparing competing traits side-by-side:
- `shortsleep_dashti` vs `shortsleep_az`
- `longsleep_dashti` vs `longsleep_az`
- `insomnia_ukb` vs `insomnia_full`

## Layout

```
config/traits.tsv                       90-trait catalog & metadata table
scripts/01_harmonize.py                 QC with an auditable ledger per trait
scripts/05_collate.py                   LDSC logs -> tidy tables + QC gate + BH-FDR + comparison summary
scripts/06_heatmap.py                   Figure 2 Heatmap (dynamic layout)
data/harmonized/*.qc.txt                every SNP dropped, and why
results/tables/h2_summary.tsv           h2 estimates, Z-scores, intercepts
results/tables/problem_comparison_summary.tsv  dual-strategy comparison table
results/tables/rg_matrix.tsv            pairwise genetic correlation matrix
results/figures/fig2_rg_heatmap.png     Figure 2 heatmap
```

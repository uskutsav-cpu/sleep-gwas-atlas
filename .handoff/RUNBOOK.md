# RUNBOOK — exact commands

All commands run from the repository root. `PY=.ldsc-env/bin/python`.

## Status
```bash
.ldsc-env/bin/python scripts/14_phase1_batch.py --status
.ldsc-env/bin/python scripts/13_build_ledger.py --summary
```

## Continue Phase 1 (resumable, safe to kill)
```bash
# one batch of 3 disease traits; repeat until the count stops rising
.ldsc-env/bin/python scripts/14_phase1_batch.py --role disease --limit 3 --min-free-gb 2.5

# uninterrupted in a human terminal:
for i in $(seq 1 40); do
  .ldsc-env/bin/python scripts/14_phase1_batch.py --role disease --limit 3 --min-free-gb 2.5
done
```

## Single trait, step by step
```bash
.ldsc-env/bin/python scripts/09_fetch_public_sources.py --trait <id>
.ldsc-env/bin/python scripts/01_harmonize.py --trait <id> --config config/traits.tsv \
    --infile data/raw/<file> --outdir data/harmonized --chunksize 1500000
bash scripts/12_phase1_run.sh <id>
```

## Rebuild the primary rg matrix and Figure 2
```bash
DIS=$(ls data/munged/*.sumstats.gz | grep -vE 'sleepdur|shortsleep|longsleep|sleepiness|chronotype' | tr '\n' ',' | sed 's/,$//')
for s in sleepdur shortsleep longsleep sleepiness chronotype; do
  .ldsc-env/bin/python ldsc/ldsc.py --rg data/munged/$s.sumstats.gz,$DIS \
      --ref-ld-chr ref/eur_w_ld_chr/ --w-ld-chr ref/eur_w_ld_chr/ --out results/logs/rg_$s
done
.ldsc-env/bin/python scripts/05_collate.py --mode rg --logdir results/logs \
    --out results/phase1/rg_primary_core_sleep.tsv
.ldsc-env/bin/python scripts/06_heatmap.py --rg results/phase1/rg_primary_core_sleep.tsv \
    --config config/traits.tsv --out figures/figure2_rg_heatmap.png
```

## Verify reproducibility of any input
```bash
.ldsc-env/bin/python scripts/09_fetch_public_sources.py --verify
```

## Audit
```bash
.ldsc-env/bin/python scripts/10_phase0_audit.py --block
.ldsc-env/bin/python scripts/10_phase0_audit.py --check-build data/raw/<file>
```

## Before Phase 2
```bash
brew install r
# then install LAVA, pin the commit, record the LD reference version
```

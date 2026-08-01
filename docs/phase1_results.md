# Phase 1 — first real LDSC results

Branch `claude/neuro-immune-cardio`. Everything below came out of LDSC run on
downloaded summary statistics; nothing here is synthetic, and every figure
carries a provenance footer tying it to the table that produced it.

---

## 1. Environment

LDSC is the maintained **CBIIT Python-3.9 fork** (`ldsc39`), not the stale
Python 2.7 original that `CLAUDE.md` warns about. It was already installed by
the other workstream, so `ldsc/`, `.ldsc-env/` and `ref/` in this repo are
**symlinks** to that shared install rather than duplicates. All three are
gitignored.

Practical consequence worth knowing: that install lives under
`~/Documents/Codex/2026-07-28/`. Deleting that directory breaks Phase 1 here.

Reference panel: `eur_w_ld_chr` (1000G Phase 3 EUR, HapMap3 SNPs), 48 files,
plus `w_hm3.snplist` as the `--merge-alleles` target. EUR panel against EUR
sumstats, per rule 3.

## 2. h² and the QC gate

Computed on the **observed scale**. Every `pop_prev` in the registry is
`UNKNOWN` by design (rule 1 — an uncited prevalence is worse than a blank), so
the liability conversion is deferred until the mentor supplies sourced values.
This does not affect the gate: Z = h²/SE is invariant to the conversion, which
rescales h² and its SE by the same constant.

| trait | h² (obs) | SE | Z | intercept | mean χ² | N_eff×h² | verdict |
|---|---|---|---|---|---|---|---|
| bmi | 0.1908 | 0.0053 | **36.0** | **1.1866** | 3.940 | 129,987 | PASS |
| mdd | 0.1938 | 0.0113 | 17.2 | 0.9962 | 1.254 | 87,182 | PASS |
| alz | 0.0119 | 0.0019 | 6.3 | 1.0302 | 1.118 | **4,551** | PASS |

Three things in that table are worth saying out loud rather than leaving in a
spreadsheet:

**BMI's intercept is 1.1866, against a gate of ≤1.20.** It passes, but only
just. That is the UK Biobank sample-overlap signature, predicted earlier from
the fact that Yengo's BMI includes UKB and so do the sleep traits. It is not
an error, but it is the number a reviewer will find, and it should be reported
rather than buried.

**Alzheimer's clears the power gate and fails the MiXeR one.** Z = 6.3 is
comfortably above 4, but N_eff × h² = 4,551 is well under CDG3's inclusion
threshold of 12,000. These are different questions — the trait is powered
enough for interpretable rg, and not powered enough to carry a MiXeR overlap
analysis. Recording both prevents the two from being conflated later.

**MDD's intercept is 0.9962**, essentially 1.0 — no detectable confounding,
which is what the public no-23andMe PGC release should look like.

## 3. Genetic correlations

18 pairs: 3 disease traits (this block) × 6 sleep/aging traits (other block,
using their munged sumstats directly). BH-FDR across all 18.

**8 survive FDR < 0.05.** Figure: `results/figures/fig2_rg_heatmap.png`.

| pair | rg | SE | FDR |
|---|---|---|---|
| mdd × frailty | **0.650** | 0.028 | <1e-6 |
| bmi × healthspan | 0.414 | 0.029 | <1e-6 |
| bmi × parental lifespan | −0.391 | 0.026 | <1e-6 |
| bmi × frailty | 0.380 | 0.017 | <1e-6 |
| mdd × healthspan | 0.300 | 0.042 | <1e-6 |
| mdd × parental lifespan | −0.270 | 0.032 | <1e-6 |
| **mdd × sleep efficiency** | **−0.171** | 0.044 | 2.3e-4 |
| **bmi × sleep efficiency** | **−0.163** | 0.025 | <1e-6 |

### The two results that are actually about sleep

Only three of the six columns are sleep phenotypes; the other three (frailty,
healthspan, parental lifespan) are aging traits from the other block. Of the
sleep columns, **sleep efficiency is the only one that correlates with
anything**: negatively with both depression (−0.17) and BMI (−0.16).

Actigraphy sleep duration and sleep midpoint/timing correlate with **nothing**
— all six of their tests are null. That is a real result and it is consistent
with the project's own prior expectation that objectively-measured duration
carries less genetic signal than the self-reported measures.

### Alzheimer's is null across the board

Not one of the six alz tests survives, and the largest is rg = −0.10 with a SE
of 0.06. Given h² = 0.0119 this is a genuine null rather than a power failure
at the rg stage, but the small h² does widen every standard error, so "no
detectable correlation" is the honest phrasing rather than "no correlation".

### Sign convention needs checking before this goes in a slide

`healthspan` correlates **positively** with BMI (+0.41) and depression (+0.30),
while `parental_lifespan` correlates **negatively** with the same two. If both
were coded in the "more/longer is healthier" direction those signs should
agree. They do not, which means at least one of the two is coded as a risk or
event rather than a duration. This is the other workstream's trait, so it is
flagged rather than changed — **confirm the coding direction before
interpreting either sign.**

## 4. Two real bugs found by running actual data

**An all-missing FRQ column silently destroyed a whole trait.** `alz` carries
no allele frequencies. Harmonize emitted `FRQ=NA` for all 8.8M SNPs, LDSC's
munge then applied `--maf-min` to that column, every value was missing, and
all 8.8M SNPs were dropped — surfacing only as the opaque error
`ValueError: No objects to concatenate`. The column is now omitted when it
holds no data, and the omission is written into the QC ledger.

**A [CDG3] filter was silently not running.** Yengo's frequency column is
named `Freq_Tested_Allele_in_HRS`, which matched no alias, so the MAF filter
did nothing on BMI while the ledger recorded "FRQ column absent". After adding
the alias it drops 155 SNPs and reports median MAF 0.215. A filter that
quietly does not run is worse than one that fails loudly.

## 5. Memory: the pipeline could not process large traits at all

`mdd` and `t2d` killed the harmonizer outright with no traceback. Cause: the
single-pass path holds the entire frame plus filtering copies in memory, and
this machine has **8 GB RAM with ~58 MB free**. A 13.5M-row file needs ~3.8 GB
for the frame alone and several times that during filtering.

Two fixes, both regression-tested to reproduce the in-memory result exactly
(BMI: 1,973,592 kept, identical either way):

- read only the columns that map to the output schema, not every column as
  Python strings;
- new `--chunksize` streaming path that filters chunk by chunk and appends.

With `--chunksize 1000000`, mdd processed its 13,554,550 rows and kept
7,149,543 (52.8%).

The chunked path makes two documented compromises, both in the ledger:
duplicate rsIDs are removed within chunks rather than across the file (LDSC's
munge removes cross-file duplicates itself), and the SNP-specific effective-N
filter uses the first chunk's median N as its reference.

## 6. Blockers

**t2d cannot enter Phase 1.** Scott 2017 curates cleanly at Phase 0 — GRCh37
verified, 12,056,347 variants — but its marker column is a bare coordinate
(`Chr:Position`) with **no rsID anywhere in the file**. LDSC's `--merge-alleles`
step needs rsIDs to merge against HapMap3, so the trait is unusable for h² or
rg without an external mapping step. This is the exact mirror of the `leptin`
blocker (rsIDs, no coordinates) and of the other block's `longevity` blocker.

The pattern is worth generalising: **Phase 0 curation and Phase 1 usability are
different bars.** A file can be public, EUR, GRCh37, genome-wide and hashed,
and still be unusable because of which identifier it happens to carry.

Route back for t2d is Mahajan 2022 through the DIAGRAM access request.

## 7. What Phase 1 still needs

- **Mentor question 2** (liability vs observed scale) before any h² table is
  presented as final. The gate results above are unaffected.
- The remaining curated traits munged. Every raw file was evicted after
  verification, so each needs re-fetching from its recorded URL — the SHA-256
  makes that exact, and it was validated in practice this session: four traits
  re-downloaded and returned byte-identical hashes.
- The sleep side is thin. Three sleep phenotypes is not an atlas; the
  self-reported traits (insomnia, chronotype, daytime sleepiness) are where the
  known signal lives and none of them is munged yet.

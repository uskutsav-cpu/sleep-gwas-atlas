# Phase 1 results — Sleep/Circadian Genetic Atlas

All values below are LDSC output on downloaded summary statistics. Nothing is
synthetic. Figures carry a provenance footer tying them to their source table.

## Scope actually completed

- **14 traits** have SNP heritability; **all 14 pass** the QC gate
  (Z >= 4 and LDSC intercept <= 1.20).
- **Primary matrix: 5 core sleep x 8 disease = 40 pairs**, BH-FDR across
  the full primary family. **15 survive FDR < 0.05.**
- h2 is on the **observed scale** — every pop_prev is UNKNOWN by design.

## SNP heritability

| trait | h2 (obs) | SE | Z | intercept | N_eff x h2 | MiXeR |
|---|---|---|---|---|---|---|
| bmi | 0.1908 | 0.0053 | 36.0 | 1.187 | 130,000 | pass |
| chronotype | 0.1208 | 0.004 | 30.2 | 1.061 | 45,630 | pass |
| pulse_pressure | 0.1024 | 0.0034 | 30.1 | 1.149 | 77,580 | pass |
| sbp | 0.1158 | 0.0039 | 29.7 | 1.194 | 87,730 | pass |
| dbp | 0.1074 | 0.0038 | 28.3 | 1.177 | 81,370 | pass |
| shortsleep | 0.0495 | 0.0021 | 23.6 | 1.027 | 15,610 | pass |
| sleepdur | 0.0664 | 0.0029 | 22.9 | 1.038 | 29,620 | pass |
| sleepiness | 0.0478 | 0.0023 | 20.8 | 1.023 | 21,610 | pass |
| mdd | 0.1938 | 0.0113 | 17.1 | 0.9962 | 87,180 | pass |
| heart_rate | 0.1468 | 0.0092 | 16.0 | 1.153 | 67,380 | pass |
| longsleep | 0.029 | 0.0021 | 13.8 | 1.015 | 3,567 | FAIL |
| crp | 0.1267 | 0.011 | 11.5 | 0.9992 | 72,920 | pass |
| alz | 0.0119 | 0.0019 | 6.3 | 1.03 | 4,551 | FAIL |
| leptin | 0.0986 | 0.016 | 6.2 | 0.9946 | 5,141 | FAIL |

## Primary genetic correlations surviving FDR

| sleep trait | disease | rg | SE | Z | FDR |
|---|---|---|---|---|---|
| longsleep | mdd | +0.404 | 0.0386 | 10.5 | 0.00e+00 |
| shortsleep | mdd | +0.340 | 0.0342 | 9.9 | 0.00e+00 |
| longsleep | leptin | +0.218 | 0.0805 | 2.7 | 1.89e-02 |
| shortsleep | bmi | +0.206 | 0.02 | 10.3 | 0.00e+00 |
| sleepiness | mdd | +0.152 | 0.0305 | 5.0 | 0.00e+00 |
| sleepiness | leptin | +0.146 | 0.0576 | 2.5 | 3.04e-02 |
| longsleep | bmi | +0.134 | 0.0246 | 5.4 | 0.00e+00 |
| sleepiness | bmi | +0.133 | 0.0201 | 6.6 | 0.00e+00 |
| longsleep | dbp | +0.104 | 0.0268 | 3.9 | 4.44e-04 |
| longsleep | heart_rate | +0.098 | 0.0305 | 3.2 | 4.73e-03 |
| sleepdur | mdd | -0.098 | 0.0315 | -3.1 | 6.33e-03 |
| longsleep | sbp | +0.096 | 0.0269 | 3.6 | 1.60e-03 |
| shortsleep | heart_rate | +0.089 | 0.0226 | 3.9 | 4.44e-04 |
| sleepdur | bmi | -0.082 | 0.0186 | -4.4 | 0.00e+00 |
| chronotype | mdd | -0.077 | 0.0256 | -3.0 | 8.00e-03 |

25 of 40 pairs did not survive.

## What the matrix shows

**The U-shape is recovered from genetics alone.** Short sleep (rg = +0.34) and
long sleep (rg = +0.40) both correlate positively with major depression, while
*continuous* sleep duration correlates negatively (rg = −0.10). That is the
classic epidemiological U-shape between sleep duration and depression,
reproduced here without any phenotypic data. The same pattern appears against
BMI: short +0.21, long +0.13, continuous −0.08.

**Morningness is protective.** Chronotype × depression is **negative**
(rg = −0.077, FDR 8.0e-3). Worth stating explicitly because the fabricated
figure this repository originally shipped showed this correlation as *positive*
(+0.14). The direction in that figure was wrong.

**Long sleep carries a cardiovascular signature.** Long sleep correlates with
diastolic BP (+0.104), systolic BP (+0.096) and resting heart rate (+0.098),
none of which appear for continuous sleep duration. Note that SBP, DBP and
pulse pressure come from one cohort of 1,028,980 and are **not independent**.

**Alzheimer's is null against every sleep trait** (max |rg| = 0.04, all FDR >
0.5). Given h² = 0.0119 the standard errors are wide, so the honest wording is
"no detectable genetic correlation at this power", not "no correlation".

## QC observations that belong in a methods discussion

**Three blood-pressure traits and BMI have intercepts of 1.15–1.19**, against a
gate of ≤ 1.20. They pass, but only just. This is the UK Biobank sample-overlap
signature and it should be reported, not buried. MDD by contrast has an
intercept of 0.9962 — the public no-23andMe PGC release behaving exactly as it
should.

**longsleep passes the power gate and fails MiXeR** (N_eff × h² = 3,567 <
12,000). Those are different questions and the tables keep them separate.
`alz` is in the same position (4,551).

## Corrections to prior project statements

`CLAUDE.md` records an expectation that short sleep and "especially long sleep"
will fail the power gate, framed as a finding rather than a bug. **Both pass**,
at Z = 23.6 and Z = 13.8. The stated prior is wrong and should be corrected.

## Not complete

This is Phase 1 for **13 of ~90 traits**, not the frozen atlas. The ledger
(`results/phase_ledger.tsv`) accounts for all 149 registry rows with an
explicit state each. Remaining work and exact commands: `.handoff/RUNBOOK.md`.
The single most valuable missing item is **insomnia**, the one core sleep
anchor without an obtainable source.

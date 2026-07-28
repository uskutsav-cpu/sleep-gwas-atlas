# provenance.md

Every field below records **where the number came from**, not just the number. Source codes:

- `MT` — main text of the paper (highest confidence)
- `TB` — a numbered table in the paper
- `FL` — figure legend or table footnote
- `MD` — methods, stated explicitly
- `IN` — inferred, not stated (reasoning given)
- `SUP` — supplementary material (not yet retrieved)
- `UNKNOWN` — not available from what we hold

Extracted 28 July 2026 from the four PDFs in `refs/pdfs/`. Nothing here is from an abstract
alone; nothing here is guessed.

---

## S02 — sleep_duration / shortsleep / longsleep

Dashti et al. 2019, *Nature Communications*. PMID 30846698. DOI 10.1038/s41467-019-08917-4.
File: `refs/pdfs/S02_Dashti2019.pdf`

| Field | Value | Source |
|---|---|---|
| n_total (continuous) | 446,118 | MT, p1 abstract + p9 methods |
| ancestry | European | MD, p9 — "453,964 subjects of European ancestry" identified by K-means on PCs; 446,118 analysed after QC |
| genome build | GRCh37 | IN — not stated. Methods (p9) say analysis used **only SNPs imputed to the HRC reference panel**; HRC v1.1 is GRCh37-only. Confident inference. |
| imputation | UK10K + 1000G Phase 3 + HRC; HRC subset used | MD, p9 |
| shortsleep ncase | **106,192** | MT, p2 |
| longsleep ncase | **34,184** | MT, p2 |
| shared ncontrol | **305,742** | MT, p2 |
| shortsleep definition | **< 7 h** | MT, p2 |
| longsleep definition | **>= 9 h** | MT, p2 |

**Correction to the earlier refs note.** These counts are in the main text, not buried in the
supplement. The supplement is not needed for this field.

---

## S03 — shortsleep / longsleep (cross-population)

Austin-Zimmerman et al. 2023, *Nature Communications*. PMID 37770476. DOI 10.1038/s41467-023-41249-y.
File: `refs/pdfs/S03_AustinZimmerman2023.pdf`

| Field | Value | Source |
|---|---|---|
| n_total | 493,142 (UK Biobank 293,037 + MVP 200,100) | TB, Table 1, p2 |
| ancestry — EUR | 445,966 (90.4%) | TB, Table 1 + MT p1 |
| ancestry — AFR | 27,785 | MT, p1 |
| ancestry — AMR | 16,250 | MT, p1 |
| ancestry — EAS | 3,141 | MT, p1 |
| genome build | GRCh37 | IN — **not stated anywhere in the paper.** Inferred from imputation panels: UKB HRC + UK10K, MVP AFR via African Genome Resources, indels via 1000G Phase 3 (MD, p9). All GRCh37. Worth confirming against the sumstats header before use. |
| shortsleep definition | **<= 5 h** | MT, p1 |
| longsleep definition | **>= 10 h** | MT, p1 |
| shortsleep ncase (cross-pop) | 62,516 | MT, p3 |
| shortsleep ncontrol (cross-pop) | 412,944 | MT, p3 |
| longsleep ncase (cross-pop) | 15,962 — **see caveat** | MT, p3 |
| shortsleep EUR | 47,054 case / 382,950 control | MT, p2 |
| shortsleep AFR | 11,352 case / 15,305 control | MT, p2 |

**Caveat on longsleep ncase.** The paper gives `n = 15,962` for the *EUR* long-sleep
meta-analysis (p2) and the same figure for the *cross-population* long-sleep meta-analysis
(p3). Table 1 gives long-sleep counts of 6,969 (UKB) + 10,713 (MVP) = 17,682 across all
ancestries, which does not reconcile with 15,962 being a cross-population total. Most likely
the cross-population long-sleep analysis is effectively EUR-only, or p3 repeats the EUR
figure. **Do not use this number until resolved** — check Supplementary Data 13/14 or the
sumstats file header. Marked `UNKNOWN` for the cross-population value until then.

---

## S11 — obstructive_sleep_apnoea

Sofer et al. 2023, *eBioMedicine*. PMID 36989840. DOI 10.1016/j.ebiom.2023.104536.
File: `refs/pdfs/S11_Sofer2023.pdf`

| Field | Value | Source |
|---|---|---|
| n_total (MVP, BMI-unadjusted) | 568,576 | MT p2 + FL, Table footnote p8 |
| n_total (MVP, BMI-adjusted) | 559,070 | FL, Table footnote p8 |
| genome build | **GRCh37 / hg37** | FL — stated twice, in the Miami plot legend (p7) and the results table footnote (p8): "genomic coordinate in genome build hg37" |
| ancestry | Multi-ancestry, stratified into White / Black / Hispanic / Asian HARE groups | MT |
| ncase / ncontrol (MVP) | `UNKNOWN` | Not in main text. Methods (p3) defer genotyping and phenotyping detail to Supplemental Notes S1/S2. Retrieve the supplement. |
| replication cohort | FinnGen v7: 27,207 cases / 280,720 controls (N = 307,927) | MT, p3 |

**Correction to the earlier refs note.** The refs spreadsheet flagged a possible GRCh38 build
because MVP data are often GRCh38. That was wrong for this paper — it states hg37 explicitly.
No liftover needed for S11.

---

## P01 — cross_disorder (CDG3)

Grotzinger et al. 2026, *Nature* 649:406–415. PMID 41372416. DOI 10.1038/s41586-025-09820-3.
File: `refs/pdfs/P01_Grotzinger2026.pdf`

| Field | Value | Source |
|---|---|---|
| n cases (total, 14 disorders) | 1,056,201 | MT, p1 abstract |
| disorders | 14 childhood- and adult-onset psychiatric disorders | MT, p1 |
| genome build | **GRCh37 / hg19** | MD, p16 — verbatim: *"Base pair location is given in genome build GRCh37/hg19 throughout the Article and its Supplementary Information."* Restated p18 for the MAGMA gene mapping. |
| ancestry | EUR-like primary; EAS-like and AFR-like in Popcorn cross-ancestry analyses | MD, p17 |
| LD reference | 1000G Phase 3, ancestry-matched | MD, pp16–18 |
| MHC handling | Excluded from all summary statistics before analysis | MD, p16 |
| per-disorder ncase/ncontrol | `UNKNOWN` | Not in main text — in Supplementary Information. Retrieve if you need per-disorder counts. |

**This resolves mentor question 1 for this dataset.** CDG3 is GRCh37/hg19, stated
unambiguously in the methods. The note in the project brief was correct.

---

## Consequences for the analysis

**The build question is now largely answered, and the answer points to GRCh37.** All four
papers examined are GRCh37 — two state it explicitly (P01, S11), two imply it unambiguously
through their imputation panels (S02, S03). Standardising on GRCh38 would mean lifting over
every one of these, which costs variants at every step and buys nothing unless a later
dataset forces it. The defensible choice is GRCh37, with the caveat that it should be
confirmed for any MVP- or TOPMed-derived dataset added later.

**The short/long sleep phenotype definitions do not match across papers, and this matters
more than the build.**

| | short | long | controls |
|---|---|---|---|
| S02 Dashti 2019 | < 7 h | >= 9 h | 7–8 h |
| S03 Austin-Zimmerman 2023 | <= 5 h | >= 10 h | 7–8 h |

These are not the same trait. Dashti's "short sleep" captures 106,192 people — roughly a
quarter of the cohort — while Austin-Zimmerman's captures a far more extreme tail. A genetic
correlation computed against one is not comparable to the same correlation computed against
the other, and pooling them or presenting them side by side without comment would be
misleading. Pick one definition per trait, state it in the methods, and if you report both,
report them as distinct phenotypes with distinct labels.

This also sharpens mentor question 2. If you go liability scale, `pop_prev` has to match the
definition you chose — the population prevalence of sleeping under 5 hours is nothing like
the prevalence of sleeping under 7 hours, and using a literature prevalence that was defined
against a different threshold silently corrupts the conversion.

---

## pop_prev provenance — unresolved for every binary trait

None of the `pop_prev` values in `config/traits.tsv` has a citation. `README.md`
records them as filled in from memory. Under rule 1 that makes all nine
unsourced, independent of whether any individual value is defensible.

A pattern was noticed in the three sleep binaries — `pop_prev` sitting very
close to the cohort case fraction — and read as evidence the values were
derived from the sample rather than from epidemiology:

| trait | sample prev | pop_prev | \|diff\| | design |
|---|---|---|---|---|
| longsleep | 0.1006 | 0.10 | 0.0006 | UKB cohort |
| t2d | 0.0858 | 0.09 | 0.0042 | ascertained case/control |
| shortsleep | 0.2578 | 0.25 | 0.0078 | UKB cohort |
| insomnia | 0.2830 | 0.30 | 0.0170 | UKB cohort |
| alz | 0.1411 | 0.05 | 0.0911 | ascertained |
| chronotype | 0.6257 | 0.50 | 0.1257 | UKB cohort |
| cad | 0.1557 | 0.06 | 0.0957 | ascertained |
| mdd | 0.3414 | 0.15 | 0.1914 | ascertained |
| scz | 0.2396 | 0.01 | 0.2296 | ascertained |

**The pattern is real but does not support that inference.** Two problems.
T2D is closer to its sample prevalence (0.0042) than two of the three sleep
traits, and 9% is a genuine literature figure for T2D — so proximity alone
does not discriminate. More fundamentally, UK Biobank is a **cohort**, not an
ascertained case/control study: for a trait dichotomised within a cohort,
sample prevalence approximates population prevalence *by construction*. Close
agreement is what these three traits should show whatever the source of the
number. The disease traits diverge because they are deliberately
case-enriched, which is a property of their design, not a sign that their
prevalences are better sourced.

So the values still need citations — but proximity is not the evidence for it,
and presenting it as such would not survive a mentor's first question.

**Blast radius is narrower than `CLAUDE.md` implies.** The observed→liability
conversion is a linear rescaling: it multiplies h² and its standard error by
the same constant. Therefore:

- **The QC gate is invariant.** Z = h²/SE is unchanged by `pop_prev`, so no
  trait passes or fails the power gate because of a wrong prevalence.
- **The rg heatmap is invariant.** rg = gcov/√(h²₁h²₂); the scaling cancels.
  Figure 2 is unaffected.
- **What *is* corrupted:** the reported liability-scale h² values in
  `h2_summary.tsv`, any cross-trait h² comparison, and the CDG3 MiXeR
  inclusion arithmetic `N_eff × h² > 12,000`.

Magnitude, for longsleep with sample prevalence held at 0.1006:

| assumed pop_prev | obs→liability multiplier |
|---|---|
| 0.03 | ×2.02 |
| 0.05 | ×2.34 |
| 0.10 | ×2.91 |
| 0.15 | ×3.31 |
| 0.20 | ×3.61 |

A 1.8-fold swing in reported h² across a plausible range of a number nobody
sourced. Not urgent enough to block Phase 1's figure; more than enough to
make an h² table indefensible.

**And the definitional problem above dominates both.** A prevalence for
"sleeping < 7 h" is not a prevalence for "sleeping ≤ 5 h". Whatever value is
chosen has to match the dichotomisation actually used, or the conversion is
wrong regardless of how well-sourced the number is.

**Mentor question 4:** where did the sleep `pop_prev` values come from, and
should these traits be on the liability scale at all?

## Still outstanding

- S03 cross-population long-sleep ncase — reconcile 15,962 against Table 1 (`UNKNOWN` for now)
- S11 MVP OSA ncase/ncontrol — in Supplemental Notes S1/S2, not retrieved
- P01 per-disorder ncase/ncontrol — in Supplementary Information, not retrieved
- Everything for the remaining 30 references — best pulled from the GWAS Catalog rather than PDFs

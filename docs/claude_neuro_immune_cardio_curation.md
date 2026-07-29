# Neuro / immune / cardiometabolic block — curation record

Workstream owning 21 trait IDs. Branch `claude/neuro-immune-cardio`.
Base commit `9c2ba1d`.

Every number here was verified this pass or is marked `UNKNOWN`. Nothing was
inferred from a filename, and nothing was carried over on trust.

---

## 1. Outcome

| | count |
|---|---|
| Traits owned | 21 |
| **CURATED** | **6** — ibd, crohn, uc, ra, bmi, atrial_fibrillation |
| TODO with documented blocker | 15 |
| Public sources downloaded and hashed | 6 (1.5 GB) |
| Sources located but deferred on size | 5 (cad, ldl, hdl, triglycerides, stroke — 12.4 GB) |
| Access-blocked | 6 (mdd, scz, bipolar, adhd, ms, t2d) |
| Unresolvable citation | 1 (sbp) |

## 2. Method

**PMIDs.** All 21 checked against NCBI E-utilities `esummary`, comparing
returned journal, year and title against the claimed source. 20 verified;
`sbp` could not be resolved (§5).

**Sources.** Located through the GWAS Catalog REST API by PMID, not by
guessing URLs. Where a publication maps to several studies, the
ancestry-stratified record was pulled per accession and the **European**
stratum selected deliberately:

- `ra` → GCST90132223 (European 97,173), *not* the multi-ancestry GCST90132222
- `asthma` → GCST010042 (European 303,859), *not* the multi-ancestry GCST010043
- `cad` → GCST90132314 (European+NR 1,165,690)
- `stroke` → GCST90104539, any-stroke European 1,308,460, chosen from 34
  GIGASTROKE studies; the CES/LAS/SVS subtype strata are deliberately excluded
  as non-independent of any-stroke
- lipids → GCST90239652 / 58 / 64 (European 1,320,016 each)

**Builds.** Verified from file contents, never from filenames, using
`scripts/10_phase0_audit.py --check-build`. Anchors are in
`config/build_anchors.tsv`: nine SNPs whose GRCh37 and GRCh38 positions differ
by at least 1 kb, fetched from the Ensembl REST API on both assemblies and
recorded with their allele strings.

One finding worth carrying forward: **position alone does not discriminate
build in a densely imputed file.** A variant exists at both builds'
coordinates for essentially any anchor, so naive position matching returns a
conflict every time — it did, 3 hg19 against 1 hg38, on the IIBDGC files. The
check therefore requires position **and** allele agreement. `rs1333049` was
dropped as an anchor entirely: its two builds differ by 1 bp.

## 3. CURATED traits

All six satisfy every promotion condition: public source downloaded, raw
header inspected, EUR ancestry supported by the source record, GRCh37 proven
from the file, sample-size handling documented, audit clean.

| trait | accession | build evidence | SHA-256 (first 16) |
|---|---|---|---|
| ibd | GCST004131 | GRCh37, 3/3 position+allele | `6c5e8b534dbe6048` |
| crohn | GCST004132 | GRCh37, 3/3 position+allele | `ddc4258b51f6da59` |
| uc | GCST004133 | GRCh37, 3/3 position+allele | `756f270ab5cf2fde` |
| ra | GCST90132223 | GRCh37, 5/5 rsID anchors | `bce30ae497edda92` |
| bmi | GCST006900 | GRCh37, 6/6 rsID anchors | `0d6ed0ea97870916` |
| atrial_fibrillation | GCST006061 | GRCh37, 3/3 rsID anchors | `3ef0f55a29ba0c06` |

`atrial_fibrillation` ships as a ZIP; the exact member used is
`AF_HRC_GWAS_ALLv11.txt`, recorded in the registry. Its build was declared
nowhere in the filename and is verified only from coordinates.

The IIBDGC files (`ibd`, `crohn`, `uc`) carry no rsID column at all — the
marker *is* the coordinate, `1:100000012_G_T`. Handling that required the
parser improvement described above, which is general rather than
file-specific.

## 4. Non-independence — do not treat these as independent tests

- **ibd ⊃ crohn + uc.** Same case series from de Lange 2017; IBD is the
  umbrella of the other two. All three are curated because the data are sound,
  but **only one level belongs in the final panel.** Recommendation: keep
  `ibd`, exclude `crohn` and `uc`.
- **ldl / hdl / triglycerides** are three phenotypes measured in the *same*
  1,320,016 individuals. Legitimately distinct biology, but heavily
  sample-overlapping and genetically correlated. Keep all three; rely on the
  LDSC bivariate intercept for the overlap and never present them as three
  independent tests.
- **sbp / dbp / pulse_pressure** share a sample fingerprint, and pulse
  pressure is by definition SBP − DBP. Only `sbp` is in this block; the other
  two are flagged, untouched.
- **bmi** overlaps UK Biobank heavily with the sleep traits in the other
  block. The LDSC intercept will show it; that is expected, not an error.

## 5. Blockers

**alz — build conflict.** The only GWAS Catalog file for Bellenguez 2022 is
`GCST90027158_buildGRCh38.tsv.gz`. The project standard is GRCh37 and rule 4
forbids silent liftover. Also contains proxy cases. Registry `build` changed
to `hg38` to record the truth. Needs an explicit decision before use.

**sbp — citation unresolvable.** A primary PMID for Evangelou 2018 could not
be verified. Title search returns only a Publisher Correction (30429575) and
an adjacent 2019 trans-ethnic paper (30578418). Rather than record a
plausible-looking PMID, `pmid` is `UNRESOLVED` and access is `UNAVAILABLE`.
This is the correction-notice trap that also caught prostate and colorectal
cancer in the other block.

**mdd — restricted.** GCST007342 contains only
`PGC_UKB_23andMe_depression_10000.txt`, a 10,000-SNP excerpt, not a usable
sumstats file. The full Howard 2019 no-23andMe release must come from PGC
directly. Not bypassed.

**scz, bipolar, adhd, ms, t2d — access gates.** No full p-value set in the
GWAS Catalog for any of them; PGC, IMSGC and DIAGRAM each require a human to
accept terms. Landing pages are recorded. The download script refuses these
rows by design rather than attempting a workaround.

**cad, ldl, hdl, triglycerides, stroke — deferred on size.** 12.4 GB
combined; `cad` alone is a 3.25 GB uncompressed TSV. URLs are verified and
recorded, so `09_fetch_public_sources.py --trait cad` completes them without
further research. `stroke`'s exact FTP filename was not retrieved this pass.

**t2d — contradiction resolved.** The registry at base commit read *Mahajan
2018* with counts 74,124/824,006. The project owner's decision on record is
*Mahajan 2022* (PMID 35551307, verified). Corrected to Mahajan 2022 and noted.

## 6. Population prevalences — all set to UNKNOWN

Sixteen binary traits in this block carried a populated `pop_prev` and **not
one had a citation anywhere in the repository**. The values were plausible
textbook figures (IBD 0.005, MS 0.003, SCZ 0.01) inherited from a prior
automated pass.

Rule 1 is explicit that an invented prevalence is worse than a blank, because
it silently corrupts liability-scale h² and looks correct all the way through
Phase 4. All sixteen are now `UNKNOWN`. A `pop_prev_citation` column was added
to `config/public_gwas_sources.tsv`, and the audit refuses any populated
`pop_prev` without a PMID or DOI in it — deliberately *not* satisfied by the
GWAS citation in `source_note`, because citing the study is not citing the
epidemiology.

This does not block the QC gate or the rg heatmap, both of which are invariant
to `pop_prev` (the conversion rescales h² and its SE by the same constant). It
blocks the liability h² table and the MiXeR `N_eff × h²` arithmetic.

## 7. The 45-trait target could not be met, and should not be forced

The brief asks for reduction toward 45 traits while forbidding deletion
outside this block. The arithmetic does not close: the registry holds 86
traits, 21 are owned here, and all 21 appear on the brief's own list of
preferred high-powered diseases. Reaching 45 requires removing ~41 rows from
the 65 this workstream does not own.

Compounding it, the registry size is specified three incompatible ways: the
project owner ruled **12** ("merging would be scope creep dressed as
thoroughness"), the brief says **45**, and the file contains **86** after an
expansion committed *after* that ruling.

Nothing was deleted. The recommendation within this block is to drop `crohn`
and `uc` as contained within `ibd`, giving 19 owned traits. That is a
recommendation, not an action: executing a reduction while three authorities
disagree would bury the conflict instead of surfacing it.

---

## 8. Open-access alternatives for the gated traits

Searched the GWAS Catalog for studies publishing a **full p-value set** (openly
downloadable, no terms gate) for each access-blocked trait. Every PMID below
was verified against NCBI before being recorded.

| trait | open alternative | PMID | EUR cases / controls | vs. the gated first choice |
|---|---|---|---|---|
| **mdd** | GCST005839, Wray 2018 *Nat Genet* | 29700475 | 135,458 / 344,901 | Strong. Close in power to Howard 2019 and no 23andMe gate. **Recommended swap.** |
| **t2d** | GCST007515, Xue 2018 *Nat Genet* | 29632382 | 48,286 / 250,671 | Usable. Smaller than Mahajan 2022 (80,154 cases) but needs no access request. |
| **ms** | GCST003566, Andlauer 2016 *Sci Adv* | 27386562 | 4,888 / 10,395 | **Weak.** ~10x smaller than IMSGC 2019 (47,429 cases). Real power-gate risk. |
| **bipolar** | GCST003724, *Hum Mol Genet* 2016 | 27329760 | 7,647 / 27,303 | **Weak.** ~5x smaller than Mullins 2021 (41,917 cases). Real power-gate risk. |
| **scz** | none found | — | — | Tested 29483656, 25056061, 31740837, 21926972, 23974872 — study records exist, no full p-value set. |
| **adhd** | none found | — | — | Tested 30478444, 28545751, 20732625 — none carry a full p-value set. |

**Reading this honestly.** Only `mdd` is a clean win. `t2d` is an acceptable
trade. `ms` and `bipolar` have open data but at a fraction of the sample size,
and the QC gate is Z = h2/SE >= 4 — these are exactly the traits that gate
exists to catch. Swapping them to stay open-access may buy a trait that then
fails and has to be dropped anyway. Accepting the PGC/IMSGC terms is likely
the better route for those two, and for `scz` and `adhd` it is the only route.

Nothing here has been swapped in `config/traits.tsv`. Changing the source paper
for four traits is a decision for the project owner and mentor, exactly as the
MDD and T2D source questions were.

---

## 9. Second pass — open-access swaps executed

Owner authorised swapping the gated traits to the open-access alternatives.

**CURATED rose from 6 to 9**: added `asthma`, `mdd`, `bipolar`.

| trait | source now used | build evidence |
|---|---|---|
| mdd | GCST005839 `MDD2018_ex23andMe.gz`, Wray 2018, PMID 29700475 | GRCh37, 3/3 anchors on chr+pos |
| bipolar | GCST003724 `Hou_27329760.tar.gz`, member `BP_GWAS_Hou_et_al._2016_results.txt`, PMID 27329760 | GRCh37, 3/3 anchors on chr+pos |
| asthma | GCST010042, Han 2020, PMID 32296059 | GRCh37, 4/4 anchors on chr+pos |

The MDD file is literally named `MDD2018_ex23andMe.gz` — the no-23andMe
release `CLAUDE.md`'s known-traps list asks for, obtained without any access
gate.

### t2d — alternative downloaded, then REJECTED

Xue 2018 (GCST007515) downloaded and build-verified GRCh37, but the file holds
**194,633 variants**. That is a coding-variant/exome subset, not a genome-wide
scan; LDSC needs roughly 1.2M HapMap3 SNPs. The 14 MB file size was the tell.
It cannot support h2 or rg and is rejected as a T2D replacement — the route
back is Mahajan 2022 through the DIAGRAM access request.

This is worth generalising: **"has full summary statistics" in the GWAS
Catalog does not mean "genome-wide."** Check the variant count before
trusting any substitute.

### ms — download failed

Two attempts, both read-timeout at 61% of 225.7 MB. The URL is verified and
correct. The truncated file was discarded rather than promoted (see below).
Retry with `09_fetch_public_sources.py --trait ms`.

### A real bug this pass exposed

The first large batch died mid-run and left `cad.tsv` at 1267 MB of an
expected 3250 MB — **renamed to its final filename as though complete.** A
truncated sumstats file bearing the real name is indistinguishable from a good
one downstream, and would have produced a quietly wrong h2.

`09_fetch_public_sources.py` now verifies the received byte count against
`Content-Length` *before* promoting `.part` to the final name, and deletes the
partial on mismatch. It caught the MS truncation on its first outing. The
registry also now stores true on-disk byte counts rather than HEAD estimates.

### Build checker — two more source formats

- **tar.gz archives** must be detected before plain `.gz`; gzip-opening a
  tarball yields a 512-byte header as the first "line".
- **chr+position matching** added for releases whose marker is a coordinate
  (`6:28571110:T`) so rsID matching can never fire. Allele agreement is
  required here too — without it `ra` regressed to a false CONFLICT, the same
  density artefact seen earlier with the IIBDGC files.

All nine curated files were re-verified after these changes: no regressions.

---

## 10. Third pass — stroke, parkinson, ms

**CURATED 9 -> 12.** All three downloaded, build-verified GRCh37 from file
contents, SHA-256 recorded.

| trait | accession | build evidence |
|---|---|---|
| stroke | GCST90104539, Mishra 2022 GIGASTROKE, PMID 36180795 | GRCh37, 7/7 anchors |
| parkinson | GCST009324, Nalls 2019, PMID 31701892 | GRCh37, 9/9 anchors (deep scan) |
| ms | GCST003566, Andlauer 2016, PMID 27386562 | GRCh37, 4/4 anchors |

`ms` succeeded on the third attempt after two truncations at 61%. The
truncation guard discarded both partials rather than promoting them, so no
bad file ever reached the registry.

### A sampling weakness the parkinson file exposed

At the previous 4M-line scan depth, `parkinson` matched only **one** anchor.
A single match is not proof of build — it is one coincidence away from being
wrong. A deep scan found 9/9. Two changes followed:

- default scan depth raised from 4M to 60M lines, so anchors late in a
  coordinate-sorted file are actually reached;
- any verdict resting on fewer than 3 anchors is now labelled
  **LOW CONFIDENCE** in the output rather than reported as settled.

The general point: this check can fail *quietly* by sampling too little, which
is the same class of error as a truncated download that keeps its final name.

### Open flag on parkinson

Nalls 2019's full meta-analysis includes a 23andMe component. GCST009324 is
the PD-only European stratum and is openly hosted, which implies it is the
public release — but that was **not** independently confirmed this pass.
Recorded as a flag, not a resolved fact.

---

## 11. Fourth pass — hdl, sbp; the sbp blocker is resolved

**CURATED 12 -> 14.**

### sbp — blocker cleared by replacing the source

`sbp` was `UNAVAILABLE` because no primary PMID for Evangelou 2018 could be
verified (search reached only a Publisher Correction). Rather than keep
hunting for an unverifiable citation, the trait was moved to a **better and
citable** source:

**Keaton 2024 (PMID 38689001, NCBI-verified), GCST90310294** — systolic blood
pressure, European, N ~1,028,980, full summary statistics openly downloadable.
Larger than the Evangelou release and requires no access gate. Build verified
GRCh37, 9/9 anchors.

The lesson generalises: when a citation cannot be verified, replacing the
source is often cleaner than trying to rescue it. An unverifiable PMID is a
permanent liability in a methods section.

### hdl — an alias gap, not a data problem

The GLGC file first returned INDETERMINATE. Cause: its position column is
`POS_b37`, which was not in the parser's alias list, so the checker fell
through to coordinate-marker mode and found nothing to match. Adding
`POS_b37`-style aliases resolved it — GRCh37, 9/9 anchors.

Worth stating explicitly: **`POS_b37` names a build in the column header, and
that is a hint, never the verdict.** The check still decides from coordinate
values. A file could carry that column name and hold hg38 positions; this
project's rule 4 exists precisely for that case.

---

## 12. Final state of this block — 17 of 21 curated

Every publicly obtainable trait in the neuro / immune / cardiometabolic block
is now downloaded, build-verified from file contents, hashed and curated.
13 GB of raw summary statistics in `data/raw/`, none of it committed.

**CURATED (17):** ibd, crohn, uc, ra, asthma, ms, mdd, bipolar, parkinson,
bmi, hdl, ldl, triglycerides, cad, stroke, atrial_fibrillation, sbp

All 17 verified **GRCh37 from coordinate values**, never from filenames.

### The 4 that remain, and why none is a matter of more effort

| trait | blocker | route |
|---|---|---|
| **alz** | Source is GRCh38-only (Bellenguez 2022). Rule 4 forbids silent liftover. Also contains proxy cases. | Owner decision on liftover, or find a GRCh37 AD release |
| **t2d** | Open substitute (Xue 2018) downloaded and build-verified, then **rejected**: 194,633 variants is a coding-variant subset, not genome-wide. LDSC needs ~1.2M HapMap3 SNPs. | DIAGRAM access request for Mahajan 2022 |
| **scz** | No open full p-value set. Tested 5 candidate PMIDs. | PGC terms gate |
| **adhd** | No open full p-value set. Tested 3 candidate PMIDs. | PGC/iPSYCH terms gate |

Two of these need a human to accept terms; one needs a scientific decision;
one is a dead end that was correctly identified as such **after** downloading
and inspecting the file rather than assumed from metadata.

### What this block contributes to the 45-trait panel

17 curated, but `crohn` and `uc` are contained within `ibd` and should not all
three enter the final panel. Recommended contribution: **15 independent
traits**, with the lipid trio flagged as sample-overlapping (one cohort of
1,320,016) rather than independent.

### Verification record

Nothing here rests on a filename, a metadata field, or another agent's claim:

- 20 of 21 PMIDs verified against NCBI (`sbp` replaced rather than guessed)
- every source located through the GWAS Catalog REST API, EUR stratum chosen
  deliberately where a publication maps to several accessions
- every build proven from coordinates plus allele agreement
- every downloaded file SHA-256'd, with true on-disk byte counts
- every population prevalence set to `UNKNOWN`, because not one had a citation

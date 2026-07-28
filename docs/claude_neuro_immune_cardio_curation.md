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

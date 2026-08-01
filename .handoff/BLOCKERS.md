# BLOCKERS

## Terminal — need a human

**insomnia (CORE SLEEP ANCHOR).** The one missing anchor, and the most
consequential blocker in the project.
- Evidence: both CTG Lab URLs for the Jansen 2019 UKB-only release
  301-redirect to `https://cncr.nl/ctg/`, returning `text/html`, 134,737 bytes
  — a landing page, not a GWAS file.
- GWAS Catalog holds 7 Jansen 2019 studies and 3 Watanabe 2022 studies, **none
  with a full p-value set**.
- The only EUR full-sumstats insomnia study is GCST90018869 (Sakaue 2021) with
  **1,402 EUR cases** — a biobank ICD-code phenotype that cannot serve as the
  core anchor.
- Needs: locate the current CTG/Jansen download, or request access.

**R is not installed.** Blocks Phase 2 (LAVA) and Phase 3 (PLACO), both R
packages. Fix: `brew install r`, then install LAVA and PLACO and pin versions.

**MATLAB / Octave not installed.** Blocks conjunctional FDR (pleioFDR requires
MATLAB). Octave may substitute but is unverified. This is the hardest Phase 3
dependency — MATLAB is commercially licensed.

## Terminal — data, not access

**t2d (Scott 2017, GCST004773).** Curates cleanly at Phase 0 — GRCh37 verified,
12,056,347 variants — but its marker column is a bare `Chr:Position` with no
rsID anywhere, so it cannot merge against HapMap3. Route back: Mahajan 2022 via
the DIAGRAM access request.

**endometriosis.** Same class: coordinate-only, no rsID column.

**epilepsy.** `standard_error` is NA for 4,138,597 of 4,875,474 variants (85%).
Harmonization correctly drops every row. Source data-quality problem — beta and
p are present but SE is absent. A human could derive SE from beta and p, but
that changes the QC provenance and should be an explicit decision.

**leptin (resolved for Phase 1, flagged for interpretation).** Has rsIDs but no
coordinate columns, so its build cannot be verified from the file. It munges
because LDSC merges on rsID, but the build claim rests on the source record
rather than on evidence in the file.

## Access gates — a human must accept terms

`scz`, `adhd`, `bipolar` upgrade, `ms` → PGC / IMSGC.
Full list with evidence: `results/phase1/exclusions_and_blockers.tsv` (24 rows).

## Operational, not scientific

Background jobs are killed when the agent's execution window closes, so bulk
processing must run as repeated short foreground batches (or be launched by a
human in a terminal, where it runs uninterrupted).

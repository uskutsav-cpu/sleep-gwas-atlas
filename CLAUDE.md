# CLAUDE.md

Persistent context for this repo. Read this before acting in any session.

## What this project is

Phase 0/1 of a Sleep/Circadian Genetic Atlas: mapping shared genetic
architecture between sleep/circadian traits and human disease. Student
research project with a mentor; work is presented in weekly meetings and
must be defensible line by line.

Methods mirror Grotzinger, Werme, Peyrot et al., "Mapping the genetic
landscape across 14 psychiatric disorders," *Nature* 649:406–415 (2026)
— referred to throughout as **CDG3**. See `methods_map.md`.

Read `README.md` and `methods_map.md` at the start of every session.

## Current state

Phase 0, not complete. No summary statistics downloaded yet. All rows in
`config/traits.tsv` are `status=TODO`. The Python scripts are tested on
synthetic data; the four shell scripts have never been executed against
real files and should be treated as first drafts.

## Non-negotiable rules

1. **Never fabricate a metadata value.** Sample sizes, case/control counts,
   and population prevalences must come from the source publication with a
   citation. `UNKNOWN` is a valid and useful answer. A plausible-looking
   invented `pop_prev` is worse than a blank, because it silently corrupts
   liability-scale h² and looks fine all the way through Phase 4.

2. **Never loosen a `[CDG3]`-tagged threshold to make a file pass.** Those
   tags mark filters taken from the Nature paper's Methods. If a file fails,
   fix the file or flag it — do not move the goalposts. Ask before changing
   any tagged threshold.

3. **Never substitute a different-ancestry LD reference panel.** The panel
   must match the ancestry of the summary statistics (EUR here). If a
   download 404s, find the correct file; don't improvise.

4. **Never silently lift over genome builds.** If a file looks like hg38
   when the config says hg19, stop and report it.

5. **Flag, don't smooth.** Anomalies — weird retention percentages, an
   intercept of 1.3, a trait with 200k SNPs when others have 1.2M — get
   surfaced explicitly, not quietly worked around.

## QC provenance

Every filter in `01_harmonize.py` is tagged `[CDG3]` (from the Nature
Methods) or `[LDSC]` (required by the tool). Current values: strand-ambiguous
removed, INFO > 0.6, MAF > 1%, MHC chr6:25–34Mb excluded, hg19, munge against
HapMap3, liability scale with sum-of-effective-N ascertainment correction.

The h² QC gate is Z = h²/SE ≥ 4 and LDSC intercept ≤ 1.2. Low Z means
underpowered; high intercept means confounding. These have different
implications and should never be reported as the same kind of failure.

## Known traps

- LDSC is Python 2.7 and its conda env file is stale. If setup fights you
  for more than ~20 minutes, switch to a maintained Python 3 fork and say
  which one.
- The Jansen 2019 insomnia GWAS includes a 23andMe component that is
  access-restricted. Use a public release that excludes it.
- `shortsleep` and especially `longsleep` are dichotomised tails of a
  continuous phenotype and are expected to fail the power gate. That is a
  finding, not a bug — report it as one.
- CDG3 used `N_eff × h² > 12,000` as the inclusion threshold for MiXeR.
  Worth computing for each sleep trait.

## Working style

- Commit after each completed step, with a message describing what changed.
- Warn before starting anything that will take more than ~30 minutes of
  compute or a multi-GB download.
- When two sources disagree on a number, report both, pick one, explain why.
- Prefer stopping and asking over guessing. This is research; a wrong number
  that looks right is the expensive failure mode, not a delay.
